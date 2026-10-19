"""Package verified source, report, data, models and representative fields."""
from pathlib import Path
import hashlib
import json
import re
import zipfile

root=Path(__file__).resolve().parent
sources=['oscillating_bl.py','test_numerics.py','requirements.txt','make_report.py',
         'budget_controls.py','adaptive_fix.py','audit_submission.py','package_submission.py',
         'README.md','REPORT.md']
files=[root/p for p in sources]
report=(root/'REPORT.md').read_text(encoding='utf-8')
for target in re.findall(r'!\[[^\]]*\]\(([^)]+)\)',report):
    assert (root/target).is_file(),f'Missing report figure: {target}'
retained=[]
for delta in ['0.100000','0.010000']:
    for method in ['baseline','vanilla','combined']:
        retained.append(f'd{delta}_s0.020_n16_seed17_{method}.npz')
retained.append('d0.010000_s0.020_n16_seed17_adaptive.npz')
for p in sorted((root/'final_results').rglob('*')):
    if not p.is_file(): continue
    if p.parent.name=='fields' and p.name not in retained: continue
    files.append(p)
for name in retained:
    assert (root/'final_results'/'fields'/name).is_file()
manifest={str(p.relative_to(root)).replace('\\','/'):{'bytes':p.stat().st_size,
          'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files}
archive=root/'Oscillating_BL_Submission.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in files: z.write(p,str(p.relative_to(root)).replace('\\','/'))
    z.writestr('MANIFEST.json',json.dumps(manifest,indent=2))
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    assert all('..' not in Path(n).parts for n in z.namelist())
print(json.dumps({'archive':str(archive),'bytes':archive.stat().st_size,'files':len(files)+1,
                  'sha256':hashlib.sha256(archive.read_bytes()).hexdigest()}))
