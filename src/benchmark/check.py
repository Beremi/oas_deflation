import sys
import os
import subprocess

executable = sys.argv[1]
cwd = os.getcwd()

#sample = subprocess.Popen([executable, os.path.join(cwd, 'master.inp')])
sample = subprocess.Popen([executable, 'master.inp'])
sample.wait()

if sample.returncode != 0:
    sys.exit(1)
else:
    sys.exit(0)
