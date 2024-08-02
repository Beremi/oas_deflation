# LD Viewer
Simple viewer for LD.out files or any text files structured as a table with a header.
- example.bat for Windows
- example ld_viewer.desktop for linux right click context menu in Dolphin
  - in KDE Plasma 5 desktops place in /home/username/.local/share/kservices5/ServiceMenus
  - in KDE Plasma 6 desktops place in /home/username/.local/share/kio/servicemenus

## Requirments:
- traits
- traitsui
- pandas
- matplotlib
- numpy

## Running from command line
- print help 
```bash 
    >>> python LD_viewer.py --help`
```
- open multiple files with preset column names to display on axis x and y
```bash
    >>> python pathto/LD_viewer.py -ld_files file1.out file2.out -y column_name_for_y -x column_name_for_x
```
- generate list of files (e.g. find LD.out files in all folders starting with results)
```bash
    >>> python pathto/LD_viewer.py -ld_files `ls results*/LD.out` -y column_name_for_y -x column_name_for_x
```

