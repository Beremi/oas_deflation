import sys
import pathlib

print(pathlib.Path(__file__))
PROJECT_DIR = pathlib.Path(__file__).resolve().parents[1]
print(PROJECT_DIR)

sys.path.append(str(PROJECT_DIR / 'src' / 'preprocessor' / 'voronoi'))
