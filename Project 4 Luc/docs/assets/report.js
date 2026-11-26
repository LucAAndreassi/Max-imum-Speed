// Written report content comes from REPORT.md; results.json supplies the 3D viewer.
import { startViewer } from "./viewer.js";

const REPORT_URL =
  "https://raw.githubusercontent.com/LucAAndreassi/Max-imum-Speed/main/Project%204%20Luc/REPORT.md";
const MARKED_URL =
  "https://cdn.jsdelivr.net/npm/marked@15.0.12/lib/marked.esm.js";
const MERMAID_URL =
  "https://cdn.jsdelivr.net/npm/mermaid@11.4.1/dist/mermaid.esm.min.mjs";

const n = (value, digits = 2) => Number(value).toFixed(digits);
const badge = pass =>
  '<span class="flag ' + (pass ? 'pass">PASS' : 'fail">FAIL') + "</span>";

async function loadReportMarkdown(report) {
  const [response, { marked }] = await Promise.all([
    fetch(REPORT_URL + "?v=" + Date.now(), { cache: "no-store" }),
    import(MARKED_URL)
  ]);
  if (!response.ok) {
    throw new Error("REPORT.md returned HTTP " + response.status + ".");
  }

  let source = (await response.text())
    .replace(/^---\r?\n[\s\S]*?\r?\n---(?:\r?\n|$)/, "")
    .replace(/<script\b[^>]*>[\s\S]*?<\/script>/gi, "");

  // Protect TeX before Markdown parses underscores, pipes, and backslashes.
  const equations = [];
  const hold = value =>
    "REPORTMATH" + (equations.push(value) - 1) + "END";
  source = source
    .replace(/\x60{3}math[^\r\n]*\r?\n([\s\S]*?)\x60{3}/g,
      (_, tex) => hold("$$\n" + tex.trim() + "\n$$"))
    .replace(/\$\$[\s\S]*?\$\$|\$[^\r\n$]+?\$/g, hold);

  const escapeHTML = value => value
    .replace(/&/g, "&amp;")
    .replace(/</g, "&lt;")
    .replace(/>/g, "&gt;");
  const html = await marked.parse(source);
  report.innerHTML = html.replace(/REPORTMATH(\d+)END/g,
    (_, i) => escapeHTML(equations[Number(i)]));

  const ids = [
    "summary", "formulation", "analysis", "verification",
    "sensitivities", "optimization", "nonlinear", "results",
    "extension", "reproducibility", "limitations"
  ];
  for (const heading of report.querySelectorAll("h2")) {
    const section = /^(\d+)\./.exec(heading.textContent.trim());
    if (section && ids[Number(section[1]) - 1]) {
      heading.id = ids[Number(section[1]) - 1];
    } else if (heading.textContent.trim().toLowerCase() === "references") {
      heading.id = "references";
    }
  }
}

function placeViewerAfterSummary(report) {
  const host = report.querySelector("#arm-viewer") ||
    document.createElement("div");
  host.id = "arm-viewer";
  host.textContent = "Loading the 3D comparison…";

  const wrapper = document.createElement("section");
  wrapper.id = "comparison";
  wrapper.className = "comparison";
  wrapper.setAttribute("aria-labelledby", "comparison-title");

  const header = document.createElement("div");
  header.className = "section-heading";
  const title = document.createElement("h2");
  title.id = "comparison-title";
  title.textContent = "Initial vs. optimized";
  const instructions = document.createElement("p");
  instructions.textContent =
    "Drag to rotate both views. Scroll or pinch to zoom. Hover or tap the tube to inspect an element.";
  header.append(title, instructions);
  wrapper.append(header, host);

  // Place the viewer after the entire summary, immediately before Section 2.
  const nextSection = report.querySelector("#formulation");
  if (nextSection) nextSection.before(wrapper);
  else report.append(wrapper);
  return host;
}

function scrollTables(report) {
  for (const table of report.querySelectorAll("table")) {
    if (table.parentElement?.classList.contains("table-scroll")) continue;
    const wrapper = document.createElement("div");
    wrapper.className = "table-scroll";
    wrapper.tabIndex = 0;
    wrapper.setAttribute("role", "region");
    wrapper.setAttribute("aria-label", "Report data table");
    table.before(wrapper);
    wrapper.append(table);
  }
}

async function enhanceReport(report) {
  // start() is called without top-level await, so this does not block page load.
  if (document.readyState !== "complete") {
    await new Promise(resolve =>
      window.addEventListener("load", resolve, { once: true }));
  }

  try {
    await window.MathJax?.startup?.promise;
    if (window.MathJax?.typesetPromise) {
      await window.MathJax.typesetPromise([report]);
    }
  } catch (error) {
    console.warn("Unable to typeset report equations:", error);
  }

  const blocks = [...report.querySelectorAll("pre > code.language-mermaid")];
  if (!blocks.length) return;
  try {
    const { default: mermaid } = await import(MERMAID_URL);
    mermaid.initialize({ startOnLoad: false, securityLevel: "strict" });
    const nodes = blocks.map(code => {
      const pre = code.parentElement;
      pre.className = "mermaid";
      pre.textContent = code.textContent;
      return pre;
    });
    await mermaid.run({ nodes });
  } catch (error) {
    console.warn("Unable to render report diagrams:", error);
  }
}

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

async function initializeViewer(host) {
  try {
    const url = new URL("results.json", import.meta.url);
    url.searchParams.set("v", Date.now());
    const response = await fetch(url, { cache: "no-store" });
    if (!response.ok) {
      throw new Error("Viewer data returned HTTP " + response.status + ".");
    }
    const comparison = formatComparison(await response.json());
    await startViewer(host, comparison, badge, n);
  } catch (error) {
    console.error("Unable to initialize the 3D viewer:", error);
    host.textContent = "Unable to load the 3D viewer. " + error.message;
  }
}

async function start() {
  const report = document.querySelector("#report");
  if (!report) throw new Error('The page needs a <main id="report"> element.');
  report.textContent = "Loading REPORT.md…";

  await loadReportMarkdown(report);
  const host = placeViewerAfterSummary(report);
  scrollTables(report);
  await Promise.all([enhanceReport(report), initializeViewer(host)]);

  if (location.hash) {
    document.getElementById(
      decodeURIComponent(location.hash.slice(1))
    )?.scrollIntoView();
  }
}

start().catch(error => {
  console.error("Unable to load the written report:", error);
  const report = document.querySelector("#report");
  if (report) report.textContent = "Unable to load REPORT.md. " + error.message;
});

