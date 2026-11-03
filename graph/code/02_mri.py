"""Export the assembled representative MRI figure without resizing its pixels."""
from pathlib import Path
import argparse
from PIL import Image

parser=argparse.ArgumentParser()
parser.add_argument('--source',type=Path,required=True,help='Assembled axial/coronal MRI PNG')
parser.add_argument('--output',type=Path,default=Path(__file__).resolve().parents[1]/'result')
args=parser.parse_args()
args.output.mkdir(parents=True,exist_ok=True)
with Image.open(args.source) as source:
    image=source.convert('RGB')
    # Keep the original pixel grid; only the export format changes.
    image.save(args.output/'Figure_2.jpg',quality=97,subsampling=0,dpi=(300,300))
