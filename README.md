# ChP–CSF biomarker analysis

R and Python code for analyses of choroid plexus volume, CSF biomarkers and cognition.

## Structure

`R` and `tools` contain shared functions; `config` contains model settings. Statistical scripts are in `statistical/code`, figure scripts in `graph/code`. The `data`, `statistical/result` and `graph/result` directories are empty in Git.

## Run locally

```sh
Rscript config/install_packages.R
python -m pip install -r requirements.txt
python run.py --data /path/analysis.csv --raw-all /path/Data_all.csv
```

This runs the statistical models and bootstrap calculations, then generates figures. Input data require authorized ADNI access. `config/local.example.json` provides the path settings. Generated outputs remain local and are ignored by Git.

SEM resampling can resume from a local cache when the input and model signatures match. Use `--no-bootstrap` only with a complete verified local cache. Use a new results directory for a different dataset.

```sh
python run.py --stage figures --data /path/analysis.csv
python run.py --stage check
```

An assembled representative MRI image can be exported with `--mri-source /path/figure.png`. No MRI image is included in this repository.
