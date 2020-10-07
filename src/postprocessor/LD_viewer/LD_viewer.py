import os
import pathlib
import argparse
os.environ["ETS_TOOLKIT"] = "qt"
os.environ["QT_API"] = "pyqt5"
from traits.api import HasStrictTraits, Instance, File, List, Enum, Button, Str, Any, Bool
from traitsui.api import Item, Group, View, HSplit, NoButtons, EnumEditor, HGroup, Handler, VGroup, Tabbed, TextEditor, Label
from traitsui.file_dialog  import open_file, TextInfo, FileInfo
from mpl_qt_editor import MPLFigureEditor
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd



class LDViewer(HasStrictTraits):
    ld_file = Any(changed=True)
    open_button = Button('Open...')
    reload_button = Button('Reload')

    figure = Instance(Figure)

    data = Any()
    labels = List([], changed=True)
    x_values = Enum(values='labels')
    y_values = Enum(values='labels')

    draw_button = Button('Draw')

    clear_axis = Bool(True)

    def __init__(self, **traits):
        super().__init__(**traits)
        if pathlib.Path('LD.out').exists():
            self.ld_file = 'LD.out'

    def _figure_default(self):
        figure = Figure()
        #figure.add_axes([0.07, 0.05, 0.85, 0.92])
        figure.add_subplot(111)
        return figure

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

    def _draw_button_fired(self):
        ax = self.figure.axes[0]
        if self.clear_axis:
            ax.cla()
        ax.plot(self.data[self.x_values], self.data[self.y_values])
        self.figure.canvas.draw()
    
    view = View(VGroup(HGroup(Item('open_button', show_label=False, id='ld_open'),
                       Item('reload_button', show_label=False),
                       Item('ld_file', id='ld_file', style='readonly')
                        ),
                HSplit(Group(
                             Item('x_values'),
                             Item('y_values'),
                             Item('clear_axis'),
                             Item('draw_button', show_label=False),
                             id='ld_group',
                             ),
                       Item('figure', show_label=False, editor=MPLFigureEditor(), dock='vertical', width=0.8),
                       show_labels=False, id='ld_hsplit', dock='tab'
                      ),
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
    parser.add_argument('ld_file',
                        #metavar='file',
                        default=None,
                        type=pathlib.Path,
                        nargs='?',
                        help='LD file')

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

    if args.ld_file:
        ld_viewer = LDViewer(ld_file=args.ld_file)
    else:
        ld_viewer = LDViewer()
    if args.x in ld_viewer.labels:
        ld_viewer.x_values = args.x
    if args.y in ld_viewer.labels:
        ld_viewer.y_values = args.y
    ld_viewer.configure_traits()
