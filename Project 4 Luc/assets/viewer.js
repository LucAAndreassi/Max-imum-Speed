import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const host=document.querySelector('#arm-viewer');
if(host) {
  try { await start(); }
  catch(e) {host.textContent='3D viewer unavailable: '+e.message+'. The report tables remain available.';}
}

async function start(){
  const response=await fetch('assets/results.json');
  if(!response.ok)throw new Error('results.json not found; run python run.py');
  const data=await response.json();
  const optimized=data.optimization?.status==='converged_preliminary_model' && data.optimization.viewer_designs;
  const cs=optimized?data.optimization.viewer_designs:data.candidates,L=data.inputs.length_m;
  const initialScale=optimized?1:8;
  host.innerHTML=`<div class="viewer-shell">
    <div class="viewer-top"><strong>${optimized?'Starting and optimized uniform tubes.':'Two sourced tubes. One load case.'}</strong><p>${optimized?'Estimated material · modeled SF ≥ 2 · shell and joints unassessed':'Preliminary candidates · no verified optimum'} · shared camera · actual hollow sections</p></div>
    <div class="viewer-controls">
      <label>Deformation <input id="deform" aria-label="Deformation scale" type="range" min="0" max="30" value="${initialScale}" step="1"><output id="scale">${initialScale}×</output></label>
      <label>State <select id="state" aria-label="Displayed state"><option value="linear">Linear bending</option><option value="nonlinear">Follower-thrust nonlinear</option></select></label>
      <button id="mode" type="button">Animate bending mode</button><button id="reset" type="button">Reset camera</button>
    </div>
    <div class="viewer-stage"><div class="viewer-labels"><span>${cs[0].name}</span><span>${cs[1].name}</span></div><div class="viewer-tooltip"></div></div>
    <div class="viewer-legend"><span>Nominal linear bending stress</span><span>0</span><span class="legend-bar"></span><span id="max-stress"></span><span>MPa · ${optimized?'estimated nominal material-strength model':'no validated failure limit'}</span></div>
    <div class="viewer-metrics"></div>
    <p class="viewer-note">Drag to orbit; scroll to zoom; hover/tap a tube to inspect. Geometry uses true tube diameters. State displacements are exaggerated by the displayed factor. Animation uses a normalized unloaded bending mode with illustrative amplitude 1% of arm length, slowed to 0.6 Hz; it does not predict vibration amplitude. The wall section remains circular; ovalization is not modeled.</p>
  </div>`;
  const stage=host.querySelector('.viewer-stage'),tip=host.querySelector('.viewer-tooltip');
  const renderer=new THREE.WebGLRenderer({antialias:true,alpha:true});
  renderer.setPixelRatio(Math.min(window.devicePixelRatio,2));
  renderer.setScissorTest(true);stage.prepend(renderer.domElement);
  const camera=new THREE.PerspectiveCamera(40,1,.001,100);
  function reset(){camera.position.set(.62,.40,.72);controls.target.set(0,.035,0);controls.update();}
  const controls=new OrbitControls(camera,renderer.domElement);controls.enableDamping=true;reset();
  const scenes=cs.map(()=>new THREE.Scene());
  const groups=scenes.map(s=>{const g=new THREE.Group();s.add(g);return g;});
  const peak=Math.max(...cs.map(c=>c.peak_bending_stress_Pa));
  host.querySelector('#max-stress').textContent=(peak/1e6).toFixed(1);
  const colorLo=new THREE.Color('#e5f3f5'),colorHi=new THREE.Color('#1b647b');
  for(let j=0;j<2;j++){
    const scene=scenes[j];scene.background=new THREE.Color('#f4f8f8');
    scene.add(new THREE.AmbientLight(0xffffff,1.2));
    const light=new THREE.DirectionalLight(0xffffff,1.6);light.position.set(1,2,2);scene.add(light);
    const support=new THREE.Mesh(new THREE.BoxGeometry(.025,.1,.1),new THREE.MeshStandardMaterial({color:0x526b77}));
    support.position.set(-L/2-.0125,0,0);scene.add(support);
    const grid=new THREE.GridHelper(1.,10,0xc6d4d9,0xe2eaec);grid.position.y=-.08;scene.add(grid);
    const arrow=new THREE.ArrowHelper(new THREE.Vector3(0,1,0),new THREE.Vector3(L/2,.02,0),.12,0xb16916,.025,.012);scene.add(arrow);
    scene.userData.arrow=arrow;
    const e=cs[j];
    const card=document.createElement('div');card.className='metric-card';
    card.innerHTML=`<h3>${e.name}</h3><dl>
      <dt>Tube mass / objective contribution</dt><dd>${(e.tube_mass_kg*1000).toFixed(2)} g</dd>
      <dt>Tube + assembly</dt><dd>${(e.total_modeled_mass_kg*1000).toFixed(2)} g</dd>
      <dt>Nominal bending stress</dt><dd>${(e.peak_bending_stress_Pa/1e6).toFixed(2)} MPa</dd>
      <dt>Linear tip deflection</dt><dd>${(e.tip_deflection_m*1000).toFixed(3)} mm</dd>
      <dt>Nonlinear FE tip (160 elements)</dt><dd>${(e.nonlinear_FEM_tip_m*1000).toFixed(3)} mm</dd>
      <dt>First bending mode</dt><dd>${e.first_bending_frequency_Hz.toFixed(2)} Hz</dd>
      <dt>Overall SF ≥ 2</dt><dd class="pending">UNRESOLVED</dd>
      ${e.optimization_summary?`<dt>Minimum modeled strength SF</dt><dd>${e.optimization_summary.minimum_modeled_safety_factor.toFixed(5)} · ${e.optimization_summary.minimum_modeled_safety_factor>=2-1e-6?'PASS':'FAIL'}</dd>`:''}
      <dt>Deflection acceptance</dt><dd class="pending">LIMIT PENDING</dd>
      <dt>Frequency acceptance</dt><dd class="pending">LIMIT PENDING</dd>
      <dt>With/without N2 optimization</dt><dd>${optimized?'COMPLETED':'PENDING'}</dd>
      </dl>`;
    host.querySelector('.viewer-metrics').append(card);
  }
  let scale=initialScale,state='linear',animate=false,lastBuild=-1;
  const pickable=[];
  function curve(c,time){
    return c.linear_shape.map((p,i)=>{
      const n=c.nonlinear_shape[i];
      const dx=state==='nonlinear'?n[0]-p[0]:0;
      const z=state==='nonlinear'?n[2]:p[2];
      const mode=animate?.01*L*c.normalized_mode[i]*Math.sin(2*Math.PI*.6*time):0;
      return new THREE.Vector3(p[0]-L/2+scale*dx,scale*(z+mode),0);
    });
  }
  function disposeGroup(g){while(g.children.length){const o=g.children[0];o.geometry?.dispose();o.material?.dispose();g.remove(o);}}
  // Two concentric surfaces plus annular ends: real hollow tube, no solid cylinder.
  function piece(points,first,last,R,ri,stress){
    const vertices=[],colors=[],indices=[],rings=last-first+1,around=32;
    for(let side=0;side<2;side++)for(let k=first;k<=last;k++){
      const tangent=points[Math.min(k+1,points.length-1)].clone().sub(points[Math.max(k-1,0)]).normalize();
      const normal=new THREE.Vector3(-tangent.y,tangent.x,0);
      const radius=side?ri:R;
      const color=colorLo.clone().lerp(colorHi,Math.min(stress[k]/peak,1));
      for(let a=0;a<=around;a++){
        const phi=a/around*2*Math.PI;
        const p=points[k].clone().addScaledVector(normal,radius*Math.cos(phi));p.z+=radius*Math.sin(phi);
        vertices.push(p.x,p.y,p.z);colors.push(color.r,color.g,color.b);
      }
    }
    const stride=around+1,offset=rings*stride;
    for(let side=0;side<2;side++)for(let k=0;k<rings-1;k++)for(let a=0;a<around;a++){
      const i=side*offset+k*stride+a,j=i+stride;
      if(side)indices.push(i,i+1,j,j,i+1,j+1);else indices.push(i,j,i+1,j,j+1,i+1);
    }
    for(const k of [0,rings-1])for(let a=0;a<around;a++){
      const i=k*stride+a,j=offset+i;
      indices.push(i,i+1,j,j,i+1,j+1);
    }
    const geo=new THREE.BufferGeometry();geo.setAttribute('position',new THREE.Float32BufferAttribute(vertices,3));geo.setAttribute('color',new THREE.Float32BufferAttribute(colors,3));geo.setIndex(indices);geo.computeVertexNormals();
    return new THREE.Mesh(geo,new THREE.MeshStandardMaterial({vertexColors:true,roughness:.55,metalness:.05,side:THREE.DoubleSide}));
  }
  function build(time){
    pickable.length=0;
    for(let j=0;j<2;j++){
      const c=cs[j],points=curve(c,time),g=groups[j];disposeGroup(g);
      for(let e=0;e<20;e++){
        const first=e*12,last=(e+1)*12;
        const mesh=piece(points,first,last,c.outer_diameter_m/2,c.outer_diameter_m/2-c.wall_m,c.sample_bending_stress_Pa);
        mesh.userData={candidate:j,element:e+1,stress:Math.max(...c.sample_bending_stress_Pa.slice(first,last+1))};
        g.add(mesh);pickable.push(mesh);
      }
      const arrow=scenes[j].userData.arrow;
      arrow.position.copy(points.at(-1));arrow.position.y+=.02;
      let angle=0;
      if(state==='nonlinear'){
        const a=c.nonlinear_shape.at(-1),b=c.nonlinear_shape.at(-2);angle=Math.atan2(a[2]-b[2],a[0]-b[0]);
      }
      arrow.setDirection(new THREE.Vector3(-Math.sin(angle),Math.cos(angle),0));
    }
  }
  host.querySelector('#deform').addEventListener('input',e=>{scale=+e.target.value;host.querySelector('#scale').textContent=scale+'×';build(performance.now()/1000);});
  host.querySelector('#state').addEventListener('change',e=>{state=e.target.value;build(performance.now()/1000);});
  host.querySelector('#mode').addEventListener('click',e=>{animate=!animate;e.target.textContent=animate?'Stop mode animation':'Animate bending mode';build(performance.now()/1000);});
  host.querySelector('#reset').addEventListener('click',reset);
  let views=[];
  const observer=new ResizeObserver(()=>{
    const w=stage.clientWidth,h=stage.clientHeight;renderer.setSize(w,h,false);
    views=w<620?[{x:0,y:h/2,w,h:h/2},{x:0,y:0,w,h:h/2}]:[{x:0,y:0,w:w/2,h},{x:w/2,y:0,w:w/2,h}];
    camera.aspect=views[0].w/views[0].h;camera.updateProjectionMatrix();
  });observer.observe(stage);
  const raycaster=new THREE.Raycaster();
  renderer.domElement.addEventListener('pointermove',event=>{
    const rect=renderer.domElement.getBoundingClientRect(),x=event.clientX-rect.left,y=rect.height-(event.clientY-rect.top);
    const j=views.findIndex(v=>x>=v.x&&x<v.x+v.w&&y>=v.y&&y<v.y+v.h);
    if(j<0)return;
    const v=views[j];raycaster.setFromCamera(new THREE.Vector2(2*(x-v.x)/v.w-1,2*(y-v.y)/v.h-1),camera);
    const hit=raycaster.intersectObjects(groups[j].children)[0];
    if(hit){const e=hit.object.userData,c=cs[j];tip.style.display='block';tip.style.left=Math.min(x+10,rect.width-245)+'px';tip.style.top=Math.max(0,rect.height-y-85)+'px';
      tip.textContent=`${c.name}\nElement ${e.element} / 20\nOD ${(c.outer_diameter_m*1000).toFixed(1)} mm · wall ${(c.wall_m*1000).toFixed(1)} mm\nNominal linear stress ≤ ${(e.stress/1e6).toFixed(2)} MPa\nFailure allowable: unresolved`;
    }else tip.style.display='none';
  });
  renderer.domElement.addEventListener('pointerleave',()=>tip.style.display='none');
  build(0);
  function render(time){
    requestAnimationFrame(render);controls.update();
    if(animate&&time-lastBuild>60){build(time/1000);lastBuild=time;}
    for(let j=0;j<views.length;j++){const v=views[j];renderer.setViewport(v.x,v.y,v.w,v.h);renderer.setScissor(v.x,v.y,v.w,v.h);renderer.render(scenes[j],camera);}
  }requestAnimationFrame(render);
}
