import { startViewer } from "./viewer.js";

const n = (value, digits = 2) => Number(value).toFixed(digits);
const badge = pass => '<span class="flag ' + (pass ? 'pass">PASS' : 'fail">FAIL') + '</span>';

function formatComparison(data) {
  const o = data.optimization;
  if (!o || o.viewer_designs?.length !== 2) {
    throw new Error("The solver export must include the starting and optimized designs.");
  }
  const study = o.study;
  const target = o.resolved_minimum_bending_frequency_Hz;
  const tolerance = o.with_nonlinear.settings.tolerances.primal;
  const designs = o.viewer_designs.map((mesh, i) => {
    const s = mesh.optimization_summary;
    const inRange = (value, bounds) => value >= bounds[0] - 1e-8 && value <= bounds[1] + 1e-8;
    const diameter = study.outer_diameter_bounds_mm;
    const wall = study.wall_bounds_mm;
    const frequency = ["First bending frequency", n(target, 4) + " Hz minimum",
      n(s.first_bending_frequency_Hz, 4) + " Hz", s.constraint_residuals.bending_frequency <= tolerance];
    const rows = [
      frequency,
      ...["linear", "nonlinear"].flatMap(state => ["tension", "compression", "shear"].map(kind => {
        const check = s.strength_checks[state][kind];
        return [(state === "linear" ? "Linear" : "N2") + " " + kind + " SF",
          "≥ " + study.required_safety_factor, n(check.safety_factor, 4), check.status === "PASS"];
      })),
      ["Outer diameter", diameter.join("–") + " mm", n(s.outer_diameter_mm, 3) + " mm", inRange(s.outer_diameter_mm, diameter)],
      ["Wall thickness", wall.map(v => n(v, 2)).join("–") + " mm", n(s.wall_mm, 3) + " mm", inRange(s.wall_mm, wall)],
      ["Inner diameter", "≥ " + study.minimum_inner_diameter_mm + " mm", n(s.inner_diameter_mm, 3) + " mm",
        s.inner_diameter_mm >= study.minimum_inner_diameter_mm - 1e-8]
    ];
    if (study.deflection_limit_mm != null) {
      rows.splice(1, 0, ["Linear tip deflection", "≤ " + study.deflection_limit_mm + " mm",
        n(s.linear_tip_deflection_mm, 4) + " mm", s.constraint_residuals.linear_deflection <= tolerance]);
    }
    return {
      mesh, summary: s, label: i ? "Optimized design" : "Starting design",
      rows, frequencyPass: frequency[3], feasible: rows.every(row => row[3]),
      mass: n(s.tube_mass_g, 3) + " g",
      allowable: s.material_properties.tensile_strength_Pa / study.required_safety_factor
    };
  });
  return {designs, target, length: data.inputs.length_m, load: o.loads, o};
}

function fillReport({designs, target, load, o}) {
  const [initial, optimum] = designs;
  const s = optimum.summary;
  const ref = o.with_nonlinear.fixed_wall_reference.point;
  const stat = (value, label) => '<div><strong>' + value + '</strong><span>' + label + '</span></div>';
  document.querySelector("#hero-stats").innerHTML =
    stat(n(s.tube_mass_g, 1) + " g", "Optimized tube mass") +
    stat(n(s.first_bending_frequency_Hz, 2) + " Hz", "First bending frequency") +
    stat(n(o.fixed_wall_mass_reduction_percent, 2) + "%", "Mass saved vs. feasible reference");
  document.querySelector("#report-summary").innerHTML =
    '<p>The ' + n(load.length_m, 4) + ' m arm carries ' + n(load.thrust_N, 2) +
    ' N combined thrust and a ' + n(load.tip_mass_kg, 3) + ' kg tip assembly at ' + load.rpm.toLocaleString("en-US") +
    ' RPM. The optimizer selects a <strong>' + n(s.outer_diameter_mm, 3) + ' mm outer diameter and ' +
    n(s.wall_mm, 3) + ' mm wall</strong>. The first bending frequency reaches the ' + n(target, 4) +
    ' Hz target, and the linear tip deflection is ' + n(s.linear_tip_deflection_mm, 3) +
    ' mm. The frequency requirement and minimum wall thickness govern the optimum.</p>' +
    '<p>The ' + initial.mass + ' starting tube reaches ' + n(initial.summary.first_bending_frequency_Hz, 3) +
    ' Hz and fails the frequency requirement. The lightest feasible fixed-wall tube on the tested grid weighs ' +
    n(ref.mass_kg * 1000, 3) + ' g at ' + n(ref.outer_diameter_mm, 1) +
    ' mm OD and ' + n(ref.wall_mm, 2) + ' mm wall. Optimizing both dimensions saves <strong>' +
    n(o.fixed_wall_mass_reduction_percent, 2) + '%</strong> against that reference. Including N2 changes the optimized mass by ' +
    n(o.nonlinear_mass_change_percent, 2) + '%.</p>';

  const responses = [
    ["Tube mass (objective)", d => d.mass],
    ["Tube and tip assembly mass", d => n(d.summary.total_modeled_mass_g, 3) + " g"],
    ["Linear tip deflection", d => n(d.summary.linear_tip_deflection_mm, 4) + " mm"],
    ["N2 tip deflection (optimization mesh)", d => n(d.summary.nonlinear_tip_deflection_mm, 4) + " mm"],
    ["Minimum nominal material SF", d => n(d.summary.minimum_modeled_safety_factor, 4)]
  ];
  const table = (name, headings, body) =>
    '<div class="table-scroll" role="region" aria-label="' + name + '" tabindex="0"><table><thead><tr>' +
    headings.map(h => '<th scope="col">' + h + '</th>').join("") + '</tr></thead><tbody>' +
    body.map(row => '<tr><th scope="row">' + row[0] + '</th>' +
      row.slice(1).map(v => '<td>' + v + '</td>').join("") + '</tr>').join("") + '</tbody></table></div>';
  const responseTable = table("Objective and responses", ["Objective and responses", initial.label, optimum.label],
    responses.map(([label, value]) => [label, ...designs.map(value)]));
  const constraintTable = table("Constraint comparison", ["Constraint", "Requirement", initial.label, optimum.label],
    initial.rows.map(([label, limit], i) => [label, limit, ...designs.map(d => d.rows[i][2] + " " + badge(d.rows[i][3]))]));
  document.querySelector("#report-results").innerHTML = responseTable + constraintTable +
    '<p>Overall structural safety remains unresolved. The nominal checks cover the estimated material model. The study leaves shell buckling, joints, laminate interaction, fatigue, and operating dynamic stability outside its scope.</p>';
}

function scrollTables() {
  for (const table of document.querySelectorAll("main > table")) {
    const wrapper = document.createElement("div");
    wrapper.className = "table-scroll";
    wrapper.tabIndex = 0;
    wrapper.setAttribute("role", "region");
    wrapper.setAttribute("aria-label", "Report data table");
    table.before(wrapper);
    wrapper.append(table);
  }
}

async function start() {
  scrollTables();
  let comparison;
  try {
    const response = await fetch(new URL("results.json", import.meta.url));
    if (!response.ok) throw new Error("Solver data returned HTTP " + response.status + ".");
    comparison = formatComparison(await response.json());
    fillReport(comparison);
  } catch (error) {
    for (const id of ["report-summary", "arm-viewer", "report-results"]) {
      document.getElementById(id).textContent = "Solver results unavailable. " + error.message;
    }
    return;
  }
  await startViewer(document.querySelector("#arm-viewer"), comparison, badge, n);
}
start();
