"""Paths shared by the plotting and bootstrap scripts."""
from pathlib import Path
import argparse

ROOT = Path(__file__).resolve().parents[1]

def arguments(needs_data=False):
    parser = argparse.ArgumentParser()
    parser.add_argument('--data', type=Path, required=needs_data)
    parser.add_argument('--results', type=Path, default=ROOT/'statistical/result')
    parser.add_argument('--output', type=Path, default=ROOT/'graph/result')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    return args

def font_path(bold=False, italic=False):
    from matplotlib.font_manager import FontProperties, findfont
    style = FontProperties(family=['Arial','DejaVu Sans'],
                           weight='bold' if bold else 'normal',
                           style='italic' if italic else 'normal')
    return findfont(style)
