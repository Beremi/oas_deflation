import sys
import os
import subprocess
import pathlib
import numpy as np

def numeric_files(directory):
    for path in sorted(directory.glob('*')):
        if not path.is_file():
            continue
        if path.name == 'version.txt':
            continue
        if path.suffix.lower() in {'.vtk', '.vtu'}:
            continue
        yield path


def load_numeric_file(path):
    rows = []
    for raw in path.read_text(encoding='utf-8', errors='replace').splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.replace('\t', ' ').split()
        try:
            rows.append([float(value) for value in parts])
        except ValueError:
            continue
    if not rows:
        return np.empty((0, 0))
    return np.asarray(rows, dtype=float)


def smoke_check_results(results_dir):
    files = list(numeric_files(results_dir))
    if not files:
        print('No numeric result files found!')
        return 1

    for path in files:
        data = load_numeric_file(path)
        if data.size and not np.all(np.isfinite(data)):
            print(f'Non-finite values found in {path.name}')
            return 1
    return 0


def simple_file_compare():
    '''Compare generated result files with references, if references exist.'''
    results_dir = pathlib.Path('results')
    check_results_dir = pathlib.Path('check_results')

    if not results_dir.exists():
        print('Results directory is missing!')
        return 1

    if not check_results_dir.exists():
        print('Reference check_results directory is missing; running finite-output smoke check.')
        return smoke_check_results(results_dir)

    files = list(numeric_files(check_results_dir))
    if not files:
        print('No numeric reference files found; running finite-output smoke check.')
        return smoke_check_results(results_dir)

    identical = []
    for f in files:
        result_path = results_dir / f.name
        if not result_path.exists():
            print(f'Result file {f.name} is missing!')
            return 1

        ld_old = load_numeric_file(f)
        ld_new = load_numeric_file(result_path)
        identical.append(
            ld_old.shape == ld_new.shape
            and np.allclose(ld_old, ld_new, rtol=1.0e-8, atol=1.0e-10)
        )
    return not all(identical)
    

if __name__ == '__main__':
    executable = sys.argv[1]
    cwd = pathlib.Path.cwd()

    #sample = subprocess.Popen([executable, os.path.join(cwd, 'master.inp')])
    sample = subprocess.Popen([executable, 'master.inp'])
    sample.wait()

    if sample.returncode != 0:
        sys.exit(1)
    else:
        simple_compare = simple_file_compare()
        print(simple_compare)
        sys.exit(simple_compare)
