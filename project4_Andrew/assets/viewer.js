import * as THREE from "three";
import { OrbitControls } from "three/addons/controls/OrbitControls.js";

const data = await (await fetch("assets/results.json")).json();
const root = document.getElementById("viewer");
const ctrl = document.getElementById("viewer-controls");
const metricsRoot = document.getElementById("viewer-metrics");
const tooltip = document.createElement("div"); tooltip.className="tooltip"; document.body.appendChild(tooltip);

root.innerHTML="";
const scene=new THREE.Scene(); scene.background=new THREE.Color(0xf8fafc);
const camera=new THREE.PerspectiveCamera(35,1,0.001,10); camera.position.set(0.48,-0.55,0.35);
const renderer=new THREE.WebGLRenderer({antialias:true}); renderer.setPixelRatio(Math.min(window.devicePixelRatio,2)); root.appendChild(renderer.domElement);
const controls=new OrbitControls(camera,renderer.domElement); controls.target.set(0.15,0,0); controls.enableDamping=true;
scene.add(new THREE.HemisphereLight(0xffffff,0x666666,2.2));
const dl=new THREE.DirectionalLight(0xffffff,1.4); dl.position.set(1,-1,2); scene.add(dl);
const grid=new THREE.GridHelper(0.8,16,0xcccccc,0xe5e7eb); grid.rotation.x=Math.PI/2; grid.position.z=-0.06; scene.add(grid);

let mode="mode1"; let colorField="mode1"; let staticScale=60; let objects=[]; let rayTargets=[]; let t0=performance.now();
const offsets={initial:new THREE.Vector3(0,-0.12,0),optimized:new THREE.Vector3(0,0.12,0)};

function seqColor(u){const x=Math.max(0,Math.min(1,u)); const c=new THREE.Color(); c.setHSL(0.60,0.72,0.88-0.52*x); return c;}
function elementModeAmplitude(design,e){
  const c=design.mode1_curves[e];
  let m=0;
  for(const d of c.displacement){m=Math.max(m,Math.hypot(d[0],d[1],d[2]));}
  return Math.min(1,m);
}
function fieldValue(design,e){return colorField==="stress"?design.stress_utilization[e]:elementModeAmplitude(design,e);}
function clearArmObjects(){for(const o of objects){scene.remove(o); if(o.geometry)o.geometry.dispose(); if(o.material)o.material.dispose();} objects=[]; rayTargets=[];}
function curvePoints(curve,amp){const pts=[]; for(let i=0;i<curve.points.length;i++){const p=curve.points[i], d=curve.displacement[i]; pts.push(new THREE.Vector3(p[0]+amp*d[0],p[1]+amp*d[1],p[2]+amp*d[2]));} return pts;}
function addArm(label,design,off,amp,curveKey){
  const curves=design[curveKey];
  curves.forEach((c,e)=>{
    const pts=curvePoints(c,amp).map(p=>p.add(off));
    const path=new THREE.CatmullRomCurve3(pts,false,"catmullrom",0.2);
    const radius=Math.max(c.radius_m,0.0015);
    const value=fieldValue(design,e);
    const geo=new THREE.TubeGeometry(path,Math.max(10,pts.length*2),radius,8,false);
    const mat=new THREE.MeshStandardMaterial({color:seqColor(value),roughness:0.65,metalness:0.08});
    const mesh=new THREE.Mesh(geo,mat);
    mesh.userData={label,element:e,radius:c.radius_m,util:design.stress_utilization[e],modeAmp:elementModeAmplitude(design,e)};
    scene.add(mesh); objects.push(mesh); rayTargets.push(mesh);
    if(colorField==="stress" && design.stress_utilization[e]>1){
      const mid=path.getPoint(0.5);
      const s=new THREE.Mesh(new THREE.SphereGeometry(radius*1.7,8,8),new THREE.MeshBasicMaterial({color:0x111111}));
      s.position.copy(mid); scene.add(s); objects.push(s);
    }
  });
  const tip=design.nodes[design.nodes.length-1];
  const motor=new THREE.Mesh(new THREE.CylinderGeometry(0.020,0.020,0.012,20),new THREE.MeshStandardMaterial({color:0x444444}));
  motor.rotation.x=Math.PI/2; motor.position.set(tip[0]+off.x,tip[1]+off.y,tip[2]+off.z+0.008); scene.add(motor); objects.push(motor);
  const rootp=design.nodes[0];
  const clamp=new THREE.Mesh(new THREE.BoxGeometry(0.018,0.05,0.05),new THREE.MeshStandardMaterial({color:0x333333}));
  clamp.position.set(rootp[0]+off.x-0.009,rootp[1]+off.y,rootp[2]+off.z); scene.add(clamp); objects.push(clamp);
  const thrust=new THREE.ArrowHelper(new THREE.Vector3(0,0,1),new THREE.Vector3(tip[0]+off.x,tip[1]+off.y,tip[2]+off.z),0.055,0x111111,0.012,0.007);
  scene.add(thrust); objects.push(thrust);
  const maneuver=new THREE.ArrowHelper(new THREE.Vector3(-1,0,0),new THREE.Vector3(tip[0]+off.x,tip[1]+off.y,tip[2]+off.z),0.035,0x555555,0.010,0.006);
  scene.add(maneuver); objects.push(maneuver);
}
function renderArms(){
  clearArmObjects();
  let amp=0, key="static_curves";
  if(mode==="static"){amp=staticScale; key="static_curves";}
  else {const phase=Math.sin((performance.now()-t0)/450); amp=0.018*phase; key=mode==="mode1"?"mode1_curves":"buckling_curves";}
  addArm("Initial",data.initial,offsets.initial,amp,key); addArm("Optimized",data.optimized,offsets.optimized,amp,key);
}
function metricCard(name,d){
  const m=d.metrics; const ok=(val,lim,type)=>type==='min'?val>=lim:val<=lim; const tag=x=>`<span class="pass">${x?'PASS':'FAIL'}</span>`;
  return `<div class="metric-card"><b>${name}</b><br>Objective - arm mass: ${(m.mass_kg*1000).toFixed(1)} g<br>True max stress: ${(m.max_stress_pa/1e6).toFixed(2)} MPa ${tag(ok(m.max_stress_pa,data.configuration.stress_limit_pa,'max'))}<br>KS stress utilization: ${m.stress_ks_util.toFixed(3)} ${tag(ok(m.stress_ks_util,1.0,'max'))}<br>Tip deflection: ${(m.tip_deflection_m*1000).toFixed(3)} mm ${tag(ok(m.tip_deflection_m,data.configuration.tip_deflection_limit_m,'max'))}<br>First frequency: ${m.frequency_hz.toFixed(2)} Hz ${tag(ok(m.frequency_hz,data.configuration.frequency_limit_hz,'min'))}<br>Buckling factor: ${m.buckling_factor.toFixed(1)} ${tag(ok(m.buckling_factor,data.configuration.buckling_factor_limit,'min'))}</div>`;
}
metricsRoot.innerHTML=metricCard("Initial",data.initial)+metricCard("Optimized",data.optimized);
ctrl.innerHTML=`<button id="staticBtn">Static deformation</button><button id="modeBtn">First mode</button><button id="buckBtn">Buckling mode</button><label>Color field <select id="field"><option value="mode1">First-mode amplitude</option><option value="stress">Stress utilization</option></select></label><label>Static exaggeration <input id="scale" type="range" min="0" max="200" value="60"></label><span id="barlabel">First-mode amplitude</span><span style="display:inline-block;width:130px;height:12px;border:1px solid #777;background:linear-gradient(to right,hsl(216,72%,88%),hsl(216,72%,36%))"></span><span id="barrange">0 to 1</span>`;
document.getElementById("staticBtn").onclick=()=>{mode="static";renderArms();};
document.getElementById("modeBtn").onclick=()=>{mode="mode1";renderArms();};
document.getElementById("buckBtn").onclick=()=>{mode="buckling";renderArms();};
document.getElementById("scale").oninput=e=>{staticScale=Number(e.target.value); if(mode==="static")renderArms();};
document.getElementById("field").onchange=e=>{colorField=e.target.value; document.getElementById("barlabel").textContent=colorField==="stress"?"Stress utilization":"First-mode amplitude"; document.getElementById("barrange").textContent=colorField==="stress"?"0 to 1 (black marker = over limit)":"0 to 1"; renderArms();};

const raycaster=new THREE.Raycaster(), mouse=new THREE.Vector2();
renderer.domElement.addEventListener("mousemove",ev=>{
  const r=renderer.domElement.getBoundingClientRect(); mouse.x=2*(ev.clientX-r.left)/r.width-1; mouse.y=1-2*(ev.clientY-r.top)/r.height;
  raycaster.setFromCamera(mouse,camera); const hit=raycaster.intersectObjects(rayTargets,false)[0];
  if(hit){const u=hit.object.userData; tooltip.style.display="block"; tooltip.style.left=(ev.clientX+12)+"px"; tooltip.style.top=(ev.clientY+12)+"px"; tooltip.textContent=`${u.label}, element ${u.element+1}: radius ${(u.radius*1000).toFixed(2)} mm, first-mode amplitude ${u.modeAmp.toFixed(3)}, stress utilization ${u.util.toFixed(3)}`;}else tooltip.style.display="none";
});
function resize(){const w=root.clientWidth,h=root.clientHeight; renderer.setSize(w,h,false); camera.aspect=w/h; camera.updateProjectionMatrix();}
window.addEventListener("resize",resize); resize(); renderArms();
let last=-1; function animate(){requestAnimationFrame(animate); if(mode!=="static"){const now=performance.now(); if(now-last>45){renderArms();last=now;}} controls.update(); renderer.render(scene,camera);} animate();
