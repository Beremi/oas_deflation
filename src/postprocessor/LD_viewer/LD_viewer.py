import os
import pathlib
import argparse
os.environ["ETS_TOOLKIT"] = "qt"
os.environ["QT_API"] = "pyqt5"
from traits.api import HasStrictTraits, Instance, File, List, Enum, Button, Str, \
                        Any, Bool, DelegatesTo, Property, Float
from traitsui.api import Item, Group, View, HSplit, NoButtons, EnumEditor, HGroup,\
                         Handler, VGroup, Tabbed, TextEditor, Label, ListEditor
from traitsui.file_dialog  import open_file, TextInfo, FileInfo, save_file
from mpl_qt_editor import MPLFigureEditor
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# HGroup Item hack
def HItem(value=None, **traits):
    return HGroup(Item(value, **traits, springy=True))

class LDFile(HasStrictTraits):
    ldfile_num = 0

    def __init__(self, **traits):
        LDFile.ldfile_num += 1
        super().__init__(**traits)

    name = Str
    ld_file = Any('', changed=True)
    open_button = Button('Open...')
    reload_button = Button('Reload')
    data = Any()
    labels = List([], changed=True)

    add_ldcurve = Button('Add new LD curve')
    ld_curves = List
    figure = Instance(Figure)

    def _add_ldcurve_fired(self):
        self.ld_curves.append(LDCurve(ldfile=self, labels=self.labels, figure=self.figure))

    def _name_default(self):
        return 'F{}'.format(self.ldfile_num)

    def _open_button_fired ( self ):
        """ Handles the user clicking the 'Open...' button.
        """
        file_name = open_file(extensions=FileInfo(), id='ld_openfile')
        if file_name != '':
            self.ld_file = file_name

    def _reload_button_fired ( self ):
        """
        """
        self.labels = self.__load_data()

    def _ld_file_changed(self):
        self.labels = self.__load_data()

    def __load_data(self):
        self.data = pd.read_csv(self.ld_file, sep='\t', header=0, index_col=False)
        return list(self.data)

    def _labels_changed(self):
        for curve in self.ld_curves:
            curve.labels = self.labels

    view = View(VGroup(HGroup(Item('open_button', show_label=False, id='ld_open'),
                       Item('reload_button', show_label=False),
                       Item('ld_file', id='ld_file'),#, style='readonly')
                        ),
                Item('name'),
                Item('add_ldcurve', show_label=False),
                Item('ld_curves', style='custom', show_label=False,
                editor=ListEditor(use_notebook=True,
                               deletable=True,
                               dock_style='tab',
                               page_name='.name')),
                ))

class LDCurve(HasStrictTraits):
    ld_num = 0

    def __init__(self, **traits):
        LDCurve.ld_num += 1
        super().__init__(**traits)

    ldfile = Instance(LDFile)
    labels = List

    name = Str

    x_values = Enum(values='labels')
    y_values = Enum(values='labels')

    x_scale = Float(1)
    y_scale = Float(1)

    figure = Instance(Figure)
    draw_this_button = Button('Draw this one')
    copy_to_clipboard = Button('Copy to clipboard')
    save_to_txt = Button('Save to txt')

    def _draw_this_button_fired(self):
        ax = self.figure.axes[0]
        self._draw_this(ax)
        self.figure.canvas.draw()

    def _draw_this(self, ax):
        #if self.clear_axis:
        #   ax.cla()
        self.draw_line(ax)
        ax.legend()

    def draw_line(self, ax):
        ax.plot(self.ldfile.data[self.x_values] * self.x_scale, self.ldfile.data[self.y_values] * self.y_scale,
                label='{}-{}'.format(self.ldfile.name, self.name), marker='o')

    def _copy_to_clipboard_fired(self):
        data = np.vstack((self.ldfile.data[self.x_values] * self.x_scale, self.ldfile.data[self.y_values] * self.y_scale)).T
        df = pd.DataFrame({self.x_values: data[:, 0], self.y_values: data[:, 1]})
        df.to_clipboard(excel=True, index=False)

    def _save_to_txt_fired(self):
        file_name = save_file(extensions=FileInfo(), id='savefile')
        if file_name == '':
            print('No filename.')
        data = np.vstack((self.ldfile.data[self.x_values] * self.x_scale, self.ldfile.data[self.y_values] * self.y_scale)).T
        df = pd.DataFrame({self.x_values: data[:, 0], self.y_values: data[:, 1]})
        df.to_csv(file_name, sep='\t', index=False)

    def _name_default(self):
        return 'LD {}'.format(self.ld_num)

    view = View(VGroup(Item('name'),
                       HGroup(HItem('x_values'),
                               Item('x_scale', show_label=False)),
                       HGroup(HItem('y_values'),
                               Item('y_scale', show_label=False)),
                       Item('draw_this_button', show_label=False),
                       Item('copy_to_clipboard', show_label=False),
                       Item('save_to_txt', show_label=False),
                       )
                )


class LDFiles(HasStrictTraits):
    figure = Instance(Figure)
    add_ldfile = Button('Add new LD file')
    ldfiles = List  #(Instance('LDFile'))

    def _add_ldfile_fired(self):
        print('XXX', self.figure)
        self.ldfiles.append(LDFile(figure=self.figure))

    view = View(
        Item('add_ldfile', show_label=False),
        Item('ldfiles', style='custom', show_label=False,
             editor=ListEditor(use_notebook=True,
                               deletable=True,
                               dock_style='tab',
                               page_name='.name')),
        resizable=True
    )


class ControlPanel(HasStrictTraits):
    ldfiles = Instance(LDFiles)
    figure = Instance(Figure)

    draw_button = Button('Draw all')
    reload_button = Button('Reload all')
    clear_button = Button('Clear')

    #clear_axis = Bool(True)

    def _ldfiles_default(self):
        print(self.figure)
        return LDFiles(figure=self.figure)

    def _clear_button_fired(self):
        ax = self.figure.axes[0]
        ax.cla()
        self.figure.canvas.draw()

    def _draw_button_fired(self):
        ax = self.figure.axes[0]
        for lds in self.ldfiles.ldfiles:
            for ld in lds.ld_curves:
                ld.draw_line(ax)
        ax.legend()
        self.figure.canvas.draw()
        print('redrawn all')

    def _reload_button_fired(self):
        for lds in self.ldfiles.ldfiles:
            lds.reload_button = True
        print('reloaded all')

    view = View((Item('@ldfiles', show_label=False),
                    HGroup(
                        Item('reload_button', show_label=False, springy=True),
                       Item('clear_button', show_label=False, springy=True),
                       Item('draw_button', show_label=False, springy=True)),),
                resizable=True
            )


class LDViewer(HasStrictTraits):

    panel = Instance(ControlPanel)
    figure = Instance(Figure)

    def _figure_default(self):
        figure = Figure()
        #figure.add_axes([0.07, 0.05, 0.85, 0.92])
        figure.add_subplot(111)
        return figure

    def _panel_default(self):
        return ControlPanel(figure=self.figure)

   # def __init__(self, **traits):
   #     super().__init__(**traits)
   #     if pathlib.Path('LD.out').exists():
   #         self.panel.ldfilesld_file = 'LD.out'

    view = View( HSplit('@panel',
                       Item('figure', show_label=False, editor=MPLFigureEditor(), dock='vertical', width=0.8),
                       show_labels=False, id='ld_hsplit', dock='tab'
                      ),
                title='LD Viewer',
                resizable=True,
                height=0.75, width=0.75,
                buttons=NoButtons,
                id='ld_viewer')


def init_parser():
    # Create the parser
    parser = argparse.ArgumentParser(description='LD Viewer',
                                     allow_abbrev=True)

    # Add the arguments
    parser.add_argument('-ld_files',
                        #metavar='file',
                        default=None,
                        type=pathlib.Path,
                        nargs='+',
                        help='LD files')

    parser.add_argument('-x',
                        type=str,
                        #nargs='?',
                        help='x axis')

    parser.add_argument('-y',
                        type=str,
                        #nargs='?',
                        help='y axis')

    # Execute the parse_args() method
    return parser.parse_args()


if __name__ == '__main__':
    args = init_parser()
    #LDFiles().configure_traits(filename='test.out')
    ld_viewer = LDViewer()
    if args.ld_files:
        for ld_idx, ld_file in enumerate(args.ld_files):
            ld_viewer.panel.ldfiles.add_ldfile = True
            ld_viewer.panel.ldfiles.ldfiles[ld_idx].ld_file = ld_file
            ld_viewer.panel.ldfiles.ldfiles[ld_idx].add_ldcurve = True
            if args.x:# in ld_viewer.labels:
                #ld_viewer.x_values = args.x
                ld_viewer.panel.ldfiles.ldfiles[ld_idx].ld_curves[0].x_values = args.x
            if args.y:# in ld_viewer.labels:
                #ld_viewer.y_values = args.y
                ld_viewer.panel.ldfiles.ldfiles[ld_idx].ld_curves[0].y_values = args.y
    ld_viewer.configure_traits()

    exit()


    args = init_parser()

    if args.ld_file:
        ld_viewer = LDViewer(ld_file=args.ld_file)
    else:
        ld_viewer = LDViewer()
    if args.x in ld_viewer.labels:
        ld_viewer.x_values = args.x
    if args.y in ld_viewer.labels:
        ld_viewer.y_values = args.y
    ld_viewer.configure_traits()
