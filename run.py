"""Run statistical analyses and figure scripts in sequence."""
from pathlib import Path
import argparse, json, os, subprocess, sys, tempfile

ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--stage',choices=['statistics','figures','all','check'],default='all')
parser.add_argument('--config',type=Path,default=ROOT/'config/local.json')
parser.add_argument('--data',type=Path)
parser.add_argument('--raw-all',type=Path)
parser.add_argument('--rscript')
parser.add_argument('--bootstrap',action=argparse.BooleanOptionalAction,default=True,help='Refit SEM bootstrap models; --no-bootstrap uses a verified local cache')
parser.add_argument('--mri-source',type=Path,help='Local assembled representative MRI image')
parser.add_argument('--results',type=Path,default=ROOT/'statistical/result')
parser.add_argument('--output',type=Path,default=ROOT/'graph/result/generated')
args=parser.parse_args()
settings=json.loads(args.config.read_text(encoding='utf8')) if args.config.exists() else {}
data=args.data or (Path(settings['analysis_csv']) if settings.get('analysis_csv') else None)
raw=args.raw_all or (Path(settings['metadata_csv']) if settings.get('metadata_csv') else None)
rscript=args.rscript or settings.get('rscript') or 'Rscript'
def call(command):
    print('Running '+Path(command[1]).name,flush=True)
    subprocess.run([str(v) for v in command],cwd=ROOT,check=True)
if args.stage=='check':
    call([sys.executable,ROOT/'tools/check_repository.py'])
    raise SystemExit()
if data is None or not data.is_file():parser.error('Supply --data or analysis_csv in config/local.json.')
data=data.resolve();args.results=args.results.resolve();args.output=args.output.resolve()
if args.stage in ['statistics','all']:
    if raw is None or not raw.is_file():parser.error('Supply --raw-all for the visit and acquisition metadata.')
    import hashlib
    signature=args.results/'sem/bootstrap_input.sha256'
    model_signature=args.results/'sem/bootstrap_signature.txt'
    digest=hashlib.sha256(data.read_bytes()).hexdigest()
    expected=[hashlib.md5(p.read_bytes()).hexdigest() for p in [data,*sorted((ROOT/'config/models').glob('*'))]]
    if not args.bootstrap:
        if not signature.exists() or signature.read_text().strip()!=digest:
            raise RuntimeError('No verified SEM bootstrap cache for this input. Rerun with --bootstrap and a new results directory.')
        if not model_signature.exists() or model_signature.read_text().splitlines()!=expected:
            raise RuntimeError('SEM bootstrap model signatures differ. Rerun with --bootstrap and a new results directory.')
    # Intermediate matched samples are kept outside the repository.
    with tempfile.TemporaryDirectory(prefix='chp_analysis_') as work:
        call([rscript,ROOT/'R/run_analysis.R','--data',data,'--raw-all',raw.resolve(),
              '--work',work,'--output',args.results,'--stage','all','--bootstrap','yes' if args.bootstrap else 'no'])
    call([sys.executable,ROOT/'statistical/code/11_bootstrap_r2.py','--data',data,'--results',args.results,'--output',args.output])
    # Cached draws must belong to the same input dataset.
    if args.bootstrap:signature.write_text(digest+'\n',encoding='ascii')
    call([sys.executable,ROOT/'statistical/code/22_summarize_sem.py','--results',args.results,'--output',args.output])
if args.stage in ['figures','all']:
    for file in ['01_flow.py','03_scatter.py']:
        call([sys.executable,ROOT/'graph/code'/file,'--data',data,'--results',args.results,'--output',args.output])
    call([sys.executable,ROOT/'graph/code/04_incremental.py','--results',args.results,'--output',args.output])
    call([sys.executable,ROOT/'graph/code/05_supplementary.py','--data',data,'--results',args.results,'--output',args.output])
    if args.mri_source:
        call([sys.executable,ROOT/'graph/code/02_mri.py','--source',args.mri_source.resolve(),'--output',args.output])
