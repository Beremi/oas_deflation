from paraview.simple import *
import glob

'''Place this file in a directory with VTU files,
edit and
run "pvpython convert_vtu_to_binary.py".
'''
def convert(fnames):
    files = glob.glob(fnames)
    files.sort()
    for f in files:
        r = XMLUnstructuredGridReader(FileName=f)
        # SaveData
        SaveData('converted_' + f, proxy=r,# Writetimestepsasfileseries=1,
            DataMode='Binary', CompressorType='ZLib', CompressionLevel='6') # LZ4, LZMA, 1-9

if __name__ == '__main__':
    convert('cells_*.vtu')
