# Publish the Project 4 report

This package is a replacement for the existing `docs/` site in Max-imum-Speed. Copy
its files into the public repository's `docs/` folder, keeping the solver in
`Project 4 Luc/`. Keep the folder structure so GitHub Pages can resolve Jekyll's
layout and the relative asset paths.

1. Commit the files to `main`. In **Settings → Pages**, choose **Deploy from a
   branch**, **main**, and **/docs**. If Pages already uses that source, the commit
   triggers a rebuild.
2. Open the Pages URL from Settings in a private window. Check the equations,
   both designs, animation, and inspection at desktop and phone widths.
3. To refresh the linked numerical results, run
   `python optimize_arm.py --update-report` from `Project 4 Luc/`, then copy the
   regenerated `assets/results.json` into `docs/assets/`. The retained report
   narrative and figures describe the supplied run. Update them if the study
   inputs or numerical results change.

The report fetches its local solver export over HTTP. Opening `index.md` in the
ordinary GitHub file view does not run the viewer. The browser needs access to
jsDelivr for the pinned Three.js and MathJax dependencies.
