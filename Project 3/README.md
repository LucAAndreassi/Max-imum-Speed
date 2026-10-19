# Oscillating boundary-layer inverse problem

Start with **REPORT.md** for the formulation, measured findings, P1--P5 diagnostics,
Path B/C experiments, recommendation memo and rubric checklist.

## Reproduce

Use Python 3.12 with the versions in requirements.txt. No GPU or PyTorch is needed.

```bash
python -m pip install -r requirements.txt
python -m unittest -v test_numerics.py
python oscillating_bl.py --suite full --seeds 17 29 41 --width 16 --collocation 160 --max-nfev 600 --threads 1 --output final_results
python budget_controls.py --results final_results
python adaptive_fix.py --results final_results --rounds 5
python audit_submission.py --results final_results
python make_report.py --results final_results
```

The main study takes a few minutes on the recorded CPU. Supplementary budget
and adaptive runs take additional time. Timing varies by machine. The archive
already contains executed results, figures, all fitted models and loss histories.
Only seven representative full field arrays are retained to keep the archive small;
the full commands regenerate the others.

The neural model uses fully trainable spatial tanh outputs with the known
sinusoidal time structure. The numerical baseline solves the finite-difference
periodic state directly and independently verifies it by BDF time integration.
These choices and their limits are explained in the report.

Before submission, the team must supply the actual interview/outreach evidence,
confirm researcher identity and sharing permission, obtain the researcher check,
send the reviewed memo, and publish/submit the repository link. These actions
are not claimed to have occurred. After inserting those details into REPORT.md,
preserve your edits before rerunning make_report.py, which regenerates the report.

To rebuild the distributable archive after all checks, run
`python package_submission.py`. All report images use relative paths so they
can render after extracting the archive or committing its contents to GitHub.
