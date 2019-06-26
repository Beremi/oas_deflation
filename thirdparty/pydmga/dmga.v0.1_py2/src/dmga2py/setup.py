from distutils.core import setup, Extension

dmga_module = Extension('dmga2py', 
                        define_macros = [('MAJOR_VERSION', '0'), ('MINOR_VERSION', '1')],
                        include_dirs = ['../', '../dmga/', '../dmga/3rd/'],
                        libraries = ['voro++'],
                        library_dirs = ['../../', '../../bin', '../dmga/3rd/voro/'],
                        sources = ['dmga2py.cpp'],
                        extra_compile_args = ['-std=c++0x', '-fPIC'])

setup (name = 'DMG Alpha Package', version = '0.1', description = 'This is C-like binding package for DMGAlpha C++ library', ext_modules = [dmga_module])
