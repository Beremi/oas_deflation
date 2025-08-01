import sys
import os
import subprocess
import pathlib
import glob
import filecmp
import numpy as np

def simple_file_compare(dir1, dir2):
    ''' Compare files in dir1 with respect to dir2
    '''
    results_dir = pathlib.Path('results')
    check_results_dir = pathlib.Path('check_results')
    
    if not (results_dir.exists() and check_results_dir.exists()):
        print('One of result directories is missing!')
        return 1
    
    files = check_results_dir.glob('*.*')
    
    identical = []
    for f in files:
        if f == 'version.txt':
            continue
        else:#if f.name.startswith('LD'):
            ld_old = np.loadtxt(f, skiprows=1, delimiter='\t')
            ld_new = np.loadtxt(results_dir / f.name, skiprows=1, delimiter='\t')
            identical.append(np.all(np.isclose(ld_old, ld_new)))
        #else:
        #    identical.append(filecmp.cmp(f, results_dir / f.name, shallow=False))
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
        simple_compare = simple_file_compare('check_results', 'results')
        print(simple_compare)
        sys.exit(simple_compare)
