Prerequisities
==============

For quick installation execute 

```bash
$ make all
```
in the command line under the main directory (where Makefile and this README file is).

Follow the instructions of the makefile script - "src/pydmga copy to your python packages" (e.g. '/home/kelidas/anaconda3/lib/python3.7/site-packages').

```bash
$ python -m site
```
(OR Add the path to src/pydmga to system PATH variable.)

Now you can test the first program:

```python
from pydmga.geometry import OrthogonalGeometry
from pydmga.container import Container
from pydmga.diagram import Diagram
geometry = OrthogonalGeometry(10, 10, 10, True, True, True);
container = Container(geometry)
container.add(2,2,2,1)
container.add(2,2,8,1)
container.add(2,8,2,2)
diagram = Diagram(container)
volume = 0.0
for cell in diagram:
    volume += cell.volume()
print("volume is", volume, "should be", 10*10*10)
```
