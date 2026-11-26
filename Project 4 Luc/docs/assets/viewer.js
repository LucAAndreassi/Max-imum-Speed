export async function startViewer(host, comparison, badge, n) {
  const { designs, length, target, o } = comparison;
  const pairs = designs.map(d =>
    '<section class="metric-card"><h3>' + d.label + '</h3><p class="dimension">' +
    (d.feasible ? badge(true) + ' All configured checks' : badge(false) + ' Configured checks') +
    '</p><dl><dt>Tube mass<small>Objective · minimize</small></dt><dd class="objective">' + d.mass + '</dd>' +
    d.rows.map(([label, limit, value, pass]) =>
      '<dt>' + label + '<small>' + limit + '</small></dt><dd>' + value + ' ' + badge(pass) + '</dd>').join("") +
    '</dl></section>').join("");
  const count = o.study.linear_elements;
  host.innerHTML =
    '<div class="viewer-tools">' +
    '<label>Displayed state<select id="shape"><option value="linear_shape">Linear bending</option>' +
    '<option value="nonlinear_shape">N2 follower thrust</option></select></label>' +
    '<label><span>Deformation <output id="scale-label" for="scale">16×</output></span>' +
    '<input id="scale" type="range" min="0" max="40" step="1" value="16"></label>' +
    '<label>Color field<select id="field"><option value="mode">First bending mode magnitude</option>' +
    '<option value="stress">Nominal bending utilization</option></select></label>' +
    '<button id="play" class="primary" type="button" aria-pressed="false">Play first mode</button>' +
    '<button id="reset" type="button">Reset view</button></div>' +
    '<div class="stage" id="stage"><div class="view-labels">' +
    designs.map(d => '<div>' + d.label + ' ' + badge(d.frequencyPass) + ' <span>' +
      n(d.summary.first_bending_frequency_Hz, 2) + ' Hz</span></div>').join("") +
    '</div></div>' +
    '<div class="legend"><div><div class="legend-bar"></div><div class="legend-ticks">' +
    '<span id="legend-low">0</span><span id="legend-high">1</span></div></div>' +
    '<div><strong id="legend-name"></strong><br><small id="legend-note"></small></div></div>' +
    '<div class="inspector"><div class="inspect-controls">' +
    '<label>Inspect design<select id="inspect-design">' +
    designs.map((d, i) => '<option value="' + i + '">' + d.label + '</option>').join("") + '</select></label>' +
    '<label>Beam station<select id="inspect-element">' +
    Array.from({length: count}, (_, i) => '<option value="' + i + '">Element ' + (i + 1) + '</option>').join("") +
    '</select></label></div><p id="element-detail"></p></div>' +
    '<div class="viewer-note"><p>Tube diameter and wall thickness are drawn at true scale. ' +
    'Deformation multiplies the displayed centerline offset. The clamp, grid, and thrust arrow are schematic. ' +
    'The rendered centerlines use solver evaluations of the cubic Hermite linear field and the integrated N2 element rotation field. ' +
    'Inspection identifies span intervals on the ' + count + '-element linear mesh.</p>' +
    '<p>Playback adds a normalized stationary first-mode displacement with an illustrative peak of 1% of span. ' +
    'Both views play at 0.6 Hz for comparison. The natural-frequency check applies to the whole arm.</p></div>' +
    '<div class="pair">' + pairs + '</div>';

  const select = id => host.querySelector("#" + id);
  const state = { scale: 16, shape: "linear_shape", field: "mode", playing: false, start: 0 };
  const inspectDesign = select("inspect-design");
  const inspectElement = select("inspect-element");
  const detail = select("element-detail");
  const stages = designs.map(d => {
    const mesh = d.mesh;
    return {
      ...d,
      stations: mesh.linear_shape.length,
      mode: mesh.normalized_mode,
      stress: mesh.sample_bending_stress_Pa.map(value => value / d.allowable)
    };
  });
  const stressMax = Math.max(...stages.flatMap(d => d.stress));
  const stressDomain = Math.max(0.01, Math.ceil(stressMax * 100) / 100);

  function updateLegend() {
    const mode = state.field === "mode";
    select("legend-name").textContent = mode ? "First bending mode · normalized magnitude" : "Nominal linear bending · tension utilization";
    select("legend-low").textContent = mode ? "0" : "0%";
    select("legend-high").textContent = mode ? "1" : n(stressDomain * 100, 0) + "%";
    select("legend-note").textContent = mode
      ? "The shape shows where the governing mode moves. The frequency target is " + n(target, 4) + " Hz for each arm."
      : "100% is the nominal tension limit at the configured safety factor. Inspection labels any over-limit interval.";
  }

  function inspect(j = +inspectDesign.value, e = +inspectElement.value) {
    const d = stages[j];
    inspectDesign.value = j;
    inspectElement.value = e;
    const first = Math.round(e * (d.stations - 1) / count);
    const last = Math.round((e + 1) * (d.stations - 1) / count);
    const util = Math.max(...d.stress.slice(first, last + 1));
    const maxMode = Math.max(...d.mode.slice(first, last + 1).map(Math.abs));
    const endDisplacement = d.mesh[state.shape][last][2];
    detail.innerHTML = '<strong>' + d.label + ', element ' + (e + 1) + '</strong> · ' +
      n(1000 * length * e / count, 1) + '–' + n(1000 * length * (e + 1) / count, 1) + ' mm from clamp · ' +
      'peak mode magnitude ' + n(maxMode, 3) + ' · end displacement ' + n(endDisplacement * 1000, 3) +
      ' mm at 1× · peak linear stress ' + n(util * d.allowable / 1e6, 2) + ' MPa · utilization ' +
      n(util * 100, 2) + '% ' + (util > 1 ? '<span class="flag fail">OVER LIMIT</span>' : '') +
      ' · global frequency ' + badge(d.frequencyPass);
  }
  updateLegend();
  inspect();
  inspectDesign.addEventListener("change", () => inspect());
  inspectElement.addEventListener("change", () => inspect());
  const stage = select("stage");
  try {
    const [THREE, { OrbitControls }] = await Promise.all([
      import("three"), import("three/addons/controls/OrbitControls.js")
    ]);
    renderComparison(THREE, OrbitControls);
  } catch (error) {
    host.querySelector(".viewer-tools").hidden = true;
    stage.classList.add("unavailable");
    stage.textContent = "3D view unavailable. Check WebGL and access to the pinned Three.js library. The element data and constraint checks are shown below.";
    console.warn("Unable to initialize 3D comparison", error);
  }

  function renderComparison(THREE, OrbitControls) {
    const renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.75));
    renderer.setScissorTest(true);
    const canvas = renderer.domElement;
    canvas.setAttribute("role", "img");
    canvas.setAttribute("aria-label", "Linked 3D views of the starting and optimized hollow tubes");
    stage.prepend(canvas);
    const camera = new THREE.PerspectiveCamera(37, 1, .002, 10);
    const controls = new OrbitControls(camera, canvas);
    controls.enablePan = false;
    controls.minDistance = length * .42;
    controls.maxDistance = length * 4;
    let views = [];
    let frame = null;
    const mobile = window.matchMedia("(max-width: 700px)");
    const low = new THREE.Color("#dceef4");
    const high = new THREE.Color("#274463");
    const color = new THREE.Color();
    const sides = 32;

    const models = stages.map(d => {
      const scene = new THREE.Scene();
      scene.background = new THREE.Color("#eaf3f5");
      scene.add(new THREE.HemisphereLight(0xffffff, 0x718694, 1.75));
      const light = new THREE.DirectionalLight(0xffffff, 1.15);
      light.position.set(.2, .7, .4);
      scene.add(light);

      const grid = new THREE.GridHelper(length * 1.44, 12, 0xb4c8d0, 0xd0dfe4);
      grid.position.y = -length * .16;
      scene.add(grid);
      const mount = new THREE.Mesh(
        new THREE.BoxGeometry(length * .045, length * .18, length * .17),
        new THREE.MeshStandardMaterial({color: 0x5f7582, roughness: .92})
      );
      mount.position.x = -length * .5225;
      scene.add(mount);
      const arrow = new THREE.ArrowHelper(new THREE.Vector3(0, 1, 0), new THREE.Vector3(), length * .18, 0xa14b1c,
        length * .035, length * .025);
      scene.add(arrow);

      const ringSize = sides;
      const offset = ringSize * d.stations;
      const positions = new Float32Array(2 * offset * 3);
      const colors = new Float32Array(positions.length);
      const indices = [];
      for (let wall = 0; wall < 2; wall++) {
        for (let k = 0; k < d.stations - 1; k++) {
          for (let a = 0; a < sides; a++) {
            const first = wall * offset + k * ringSize;
            const i = first + a;
            const j = first + (a + 1) % sides;
            if (wall) indices.push(i, i + ringSize, j, j, i + ringSize, j + ringSize);
            else indices.push(i, j, i + ringSize, j, j + ringSize, i + ringSize);
          }
        }
      }
      for (const k of [0, d.stations - 1]) {
        for (let a = 0; a < sides; a++) {
          const i = k * ringSize + a;
          const j = k * ringSize + (a + 1) % sides;
          if (k === 0) indices.push(i, i + offset, j, j, i + offset, j + offset);
          else indices.push(i, j, i + offset, j, j + offset, i + offset);
        }
      }
      const geometry = new THREE.BufferGeometry();
      geometry.setIndex(indices);
      geometry.setAttribute("position", new THREE.BufferAttribute(positions, 3).setUsage(THREE.DynamicDrawUsage));
      geometry.setAttribute("color", new THREE.BufferAttribute(colors, 3));
      const tube = new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({vertexColors: true, roughness: .69}));
      scene.add(tube);
      return {d, scene, arrow, tube, geometry, offset, points: new Float64Array(d.stations * 2)};
    });

    function updateColors() {
      for (const { d, geometry, offset } of models) {
        const attr = geometry.attributes.color;
        for (let k = 0; k < d.stations; k++) {
          const value = state.field === "mode" ? Math.abs(d.mode[k]) : d.stress[k] / stressDomain;
          color.copy(low).lerp(high, Math.min(1, value));
          for (let a = 0; a < sides; a++) {
            const i = k * sides + a;
            attr.setXYZ(i, color.r, color.g, color.b);
            attr.setXYZ(i + offset, color.r, color.g, color.b);
          }
        }
        attr.needsUpdate = true;
      }
    }

    function updateShape(modeOffset) {
      for (const {d, geometry, offset, points, arrow} of models) {
        const raw = d.mesh[state.shape];
        for (let k = 0; k < d.stations; k++) {
          const base = d.mesh.linear_shape[k][0];
          points[2 * k] = base - length / 2 + state.scale * (raw[k][0] - base);
          points[2 * k + 1] = state.scale * raw[k][2] + modeOffset * d.mode[k];
        }
        const attr = geometry.attributes.position;
        const outer = d.mesh.outer_diameter_m / 2;
        for (let k = 0; k < d.stations; k++) {
          const a = 2 * Math.max(0, k - 1);
          const b = 2 * Math.min(d.stations - 1, k + 1);
          const theta = Math.atan2(points[b + 1] - points[a + 1], points[b] - points[a]);
          const nx = -Math.sin(theta), ny = Math.cos(theta);
          for (let j = 0; j < sides; j++) {
            const phi = j * 2 * Math.PI / sides;
            const c = Math.cos(phi), s = Math.sin(phi);
            const i = k * sides + j;
            const x = points[2 * k], y = points[2 * k + 1];
            const inner = outer - d.mesh.wall_m;
            attr.setXYZ(i, x + outer * nx * c, y + outer * ny * c, outer * s);
            attr.setXYZ(i + offset, x + inner * nx * c, y + inner * ny * c, inner * s);
          }
        }
        attr.needsUpdate = true;
        geometry.computeVertexNormals();
        geometry.computeBoundingSphere();
        const k = d.stations - 1;
        const theta = state.shape === "linear_shape" ? 0 :
          Math.atan2(raw[k][2] - raw[k - 1][2], raw[k][0] - raw[k - 1][0]);
        arrow.position.set(points[2 * k], points[2 * k + 1], 0);
        arrow.setDirection(new THREE.Vector3(-Math.sin(theta), Math.cos(theta), 0));
      }
    }

    function requestDraw() {
      if (frame == null) frame = requestAnimationFrame(draw);
    }
    function useView(view) {
      camera.aspect = view.width / view.height;
      camera.updateProjectionMatrix();
    }
    function draw(time) {
      frame = null;
      if (state.playing) {
        updateShape(length * .01 * Math.sin((time - state.start) / 1000 * 2 * Math.PI * .6));
      }
      views.forEach((view, i) => {
        useView(view);
        renderer.setViewport(view.x, view.y, view.width, view.height);
        renderer.setScissor(view.x, view.y, view.width, view.height);
        renderer.render(models[i].scene, camera);
      });
      if (state.playing) requestDraw();
    }
    function reset() {
      const factor = mobile.matches ? 1.13 : 1;
      camera.position.set(length * .77 * factor, length * .53 * factor, length * 1.33 * factor);
      controls.target.set(0, length * .012, 0);
      controls.update();
      requestDraw();
    }
    function resize() {
      const width = stage.clientWidth, height = stage.clientHeight;
      if (!width || !height) return;
      renderer.setSize(width, height, false);
      const w = Math.floor(width / 2), h = Math.floor(height / 2);
      views = mobile.matches
        ? [{x: 0, y: h, width, height: height - h}, {x: 0, y: 0, width, height: h}]
        : [{x: 0, y: 0, width: w, height}, {x: w, y: 0, width: width - w, height}];
      requestDraw();
    }
    controls.addEventListener("change", requestDraw);
    select("scale").addEventListener("input", event => {
      state.scale = +event.target.value;
      select("scale-label").textContent = state.scale + "×";
      if (!state.playing) updateShape(0);
      requestDraw();
    });
    select("shape").addEventListener("change", event => {
      state.shape = event.target.value;
      if (!state.playing) updateShape(0);
      inspect();
      requestDraw();
    });
    select("field").addEventListener("change", event => {
      state.field = event.target.value;
      updateColors();
      updateLegend();
      requestDraw();
    });
    select("play").addEventListener("click", event => {
      state.playing = !state.playing;
      state.start = performance.now();
      event.target.setAttribute("aria-pressed", state.playing);
      event.target.textContent = state.playing ? "Pause first mode" : "Play first mode";
      if (!state.playing) updateShape(0);
      requestDraw();
    });
    select("reset").addEventListener("click", reset);

    const ray = new THREE.Raycaster();
    function pick(event) {
      const rect = canvas.getBoundingClientRect();
      const x = event.clientX - rect.left, y = rect.bottom - event.clientY;
      const j = views.findIndex(v => x >= v.x && x < v.x + v.width && y >= v.y && y < v.y + v.height);
      if (j < 0) return;
      const view = views[j];
      useView(view);
      ray.setFromCamera(new THREE.Vector2(2 * (x - view.x) / view.width - 1, 2 * (y - view.y) / view.height - 1), camera);
      const m = models[j];
      const hit = ray.intersectObject(m.tube)[0];
      if (!hit) return;
      const bands = m.d.stations - 1;
      const bandFaces = bands * sides * 2;
      const face = hit.faceIndex;
      const k = face < 2 * bandFaces ? Math.floor((face % bandFaces) / (sides * 2)) :
        face < 2 * bandFaces + sides * 2 ? 0 : bands - 1;
      const e = Math.min(count - 1, Math.floor(k * count / bands));
      if (+inspectDesign.value !== j || +inspectElement.value !== e) inspect(j, e);
    }
    canvas.addEventListener("pointermove", e => {
      if (e.pointerType === "mouse" && e.buttons === 0) pick(e);
    });
    canvas.addEventListener("click", pick);
    new ResizeObserver(resize).observe(stage);
    mobile.addEventListener("change", resize);
    updateColors();
    updateShape(0);
    reset();
    resize();
  }
}
