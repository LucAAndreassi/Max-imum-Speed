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
  async function start() {
  await loadReportMarkdown();
  scrollTables();

  const host = document.querySelector("#arm-viewer");
  if (!host) return;

  try {
    const response = await fetch(
      new URL("results.json", import.meta.url)
    );

    if (!response.ok) {
      throw new Error("Viewer data returned HTTP " + response.status);
    }

    const comparison = formatComparison(await response.json());
    await startViewer(host, comparison, badge, n);
  } catch (error) {
    host.textContent = "Unable to load the 3D viewer. " + error.message;
  }
}

start().catch(error => {
  document.querySelector("#report").textContent =
    "Unable to load REPORT.md. " + error.message;
});
