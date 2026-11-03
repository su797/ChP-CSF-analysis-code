"""Check that tracked files contain code and configuration only."""
from pathlib import Path
import ast,csv,subprocess
root=Path(__file__).resolve().parents[1]
for p in root.rglob('*.py'):
    if '.git' not in p.parts:ast.parse(p.read_text(encoding='utf8'),filename=str(p))
files=subprocess.check_output(['git','ls-files'],cwd=root,text=True).splitlines()
for name in files:
    p=Path(name)
    if any(part in {'data','result','results','aggregate_results','bootstrap_replicates'} for part in p.parts):
        assert p.name=='.gitkeep',name
    assert p.suffix.lower() not in {'.png','.jpg','.jpeg','.svg','.pdf','.pptx','.docx','.zip','.npy','.npz','.rds','.xlsx','.html'},name
    if p.suffix=='.csv':assert name=='config/transformation_plan.csv',name
with (root/'config/transformation_plan.csv').open(encoding='utf8') as f:
    assert next(csv.reader(f))==['variable','analysis_var','transform_method']
print('Passed: tracked files contain code and configuration only.')
