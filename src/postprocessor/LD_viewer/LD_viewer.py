import os
import time
import pathlib
import argparse
import logging
import datetime
os.environ["ETS_TOOLKIT"] = "qt"
os.environ["QT_API"] = "pyqt5"
from traits.api import HasStrictTraits, Instance, File, List, Enum, Button, Str, \
                        Any, Bool, DelegatesTo, Property, Float, Tuple
from traitsui.api import Item, Group, View, HSplit, NoButtons, EnumEditor, HGroup,\
                         Handler, VGroup, Tabbed, TextEditor, Label, ListEditor, TreeEditor, TreeNode
from traitsui.menu \
    import Menu, Action, Separator
from traitsui.qt.tree_editor \
    import NewAction, CopyAction, CutAction, \
           PasteAction, DeleteAction, RenameAction
from traitsui.message import Message
from traitsui.file_dialog  import open_file, TextInfo, FileInfo, save_file
from mpl_qt_editor import MPLFigureEditor
from matplotlib.figure import Figure
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import mplcursors
import asyncio
from threading import Thread, Event

class TimeFilter(logging.Filter):

    def filter(self, record):
        try:
          last = self.last
        except AttributeError:
          last = record.relativeCreated

        delta = datetime.datetime.fromtimestamp(record.relativeCreated/1000.0) - datetime.datetime.fromtimestamp(last/1000.0)

        record.relative = '{0:.6f}'.format(delta.seconds + delta.microseconds/1000000.0)

        self.last = record.relativeCreated
        return True

logging.basicConfig(level=logging.INFO)
fmt = logging.Formatter(fmt="%(asctime)s %(levelname)s (%(relative)ss) - %(message)s")
log = logging.getLogger()
[hndl.addFilter(TimeFilter()) for hndl in log.handlers]
[hndl.setFormatter(fmt) for hndl in log.handlers]

# HGroup Item hack
def HItem(value=None, **traits):
    return HGroup(Item(value, **traits, springy=True))


def error(message="", title="Message", buttons=["OK", "Cancel"], parent=None):
    """ Displays a message to the user as a modal window with the specified
    title and buttons.

    If *buttons* is not specified, **OK** and **Cancel** buttons are used,
    which is appropriate for confirmations, where the user must decide whether
    to proceed. Be sure to word the message so that it is clear that clicking
    **OK** continues the operation.
    """
    msg = Message(message=message)
    ui = msg.edit_traits(
        parent=parent,
        view=View(
            ["message~", "|<>"], title=title, buttons=buttons, kind="modal",
            resizable=True,
            #height=150,
            #width=300,
            scrollable=True
        ),
    )
    return ui.result


class LDFile(HasStrictTraits):
    ldfile_num = 0

    def __init__(self, **traits):
        LDFile.ldfile_num += 1
        super().__init__(**traits)

    name = Str
    ld_file = Str('', changed=True, enter_set=True, auto_set=False)
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
            self.ld_file = str(pathlib.Path(file_name).absolute())

    def _reload_button_fired ( self ):
        """
        """
        self.labels = self.__load_data()

    def _ld_file_changed(self):
        self.labels = self.__load_data()

    def __load_data(self):
        try:
            logging.debug(f'Loading: {self.ld_file}')
            self.data = pd.read_csv(self.ld_file, sep='\t', header=0, index_col=False)
        except BaseException as err:
            # https://stackoverflow.com/questions/1483429/how-to-print-an-exception-in-python
            message = 'The file {} cannot be loaded\n'.format(self.ld_file)
            message += str(err)
            error(message, 'Error')
            self.data = np.array([])
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
        # if self.figure_settings.show_legend:
        #     ax.legend()
        # if self.figure_settings.set_xlim:
        #     ax.set_xlim(self.figure_settings.x_limits)
        # if self.figure_settings.set_ylim:
        #     ax.set_ylim(self.figure_settings.y_limits)

    def draw_line(self, ax):
        x_to_draw = self.ldfile.data[self.x_values]
        y_to_draw = self.ldfile.data[self.y_values]
        # if self.figure_settings.ignore_last_step:
        #     x_to_draw = x_to_draw.iloc[:-1]
        #     y_to_draw = y_to_draw.iloc[:-1]
        return ax.plot(x_to_draw * self.x_scale, y_to_draw * self.y_scale,
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
    name = 'ldfiles'
    figure = Instance(Figure)
    add_ldfile = Button('Add new LD file')
    ldfiles = List(LDFile) #(Instance('LDFile'))

    def _add_ldfile_fired(self):
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
             

class FigureSettings(HasStrictTraits):
    ignore_last_step = Bool(False)
    show_legend = Bool()
    set_xlim = Bool()
    x_limits = Tuple(np.nan, np.nan)
    set_ylim = Bool()
    y_limits = Tuple(np.nan, np.nan)
    view = View(
        Item('ignore_last_step', show_label=True),
        Item('show_legend', show_label=True),
        HGroup(Item('set_xlim', show_label=False),
               Item('x_limits', show_label=True)),
        HGroup(Item('set_ylim', show_label=False),
               Item('y_limits', show_label=True)),
        resizable=True
    )
        

class ReloadingThread(Thread):

    def run(self):
        while not self.event.is_set():
            self.controler._reload_button_fired()
            self.controler._clear_button_fired()
            self.controler._draw_button_fired()
            self.event.wait(timeout=self.controler.time_to_sleep)
            #time.sleep(self.controler.time_to_sleep)
        return



ldfile_view = View(VGroup(HGroup(Item('open_button', show_label=False, id='ld_open'),
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
                
# Tree editor
tree_editor = TreeEditor(
    nodes = [
        TreeNode( node_for  = [ LDFiles ],
                  auto_open = True,
                  children = 'ldfiles',
                  label     = '=LD Files',
                  add = [ LDFile ],
                  # menu=Menu( NewAction,
                  #            Separator(),
                  #            #def_title_action,
                  #            #dept_action,
                  #            Separator(),
                  #            CopyAction,
                  #            CutAction,
                  #            PasteAction,
                  #            Separator(),
                  #            DeleteAction,
                  #            Separator(),
                  #            RenameAction),
                  view =  View()),
        TreeNode( node_for  = [ LDFile ],
                  auto_open = True,
                  label     = 'name',
                  menu=Menu( NewAction,
                             Separator(),
                             #def_title_action,
                             #dept_action,
                             Separator(),
                             CopyAction,
                             CutAction,
                             PasteAction,
                             Separator(),
                             DeleteAction,
                             Separator(),
                             RenameAction),
                  view =  ldfile_view)
    ],
    #hide_root=True,
    orientation="vertical",
)
                  

class ControlPanel(HasStrictTraits):
    ldfiles = Instance(LDFiles)
    figure_settings = Instance(FigureSettings, ())
    figure = Instance(Figure)

    draw_button = Button('Draw all')
    reload_button = Button('Reload all')
    clear_button = Button('Clear')
    
    get_current_limits = Button('Get current limits')
        
    reloading_thread = Instance(Thread)
    event = Instance(Event, ())
    
    time_to_sleep = Float(60)
    start_stop_reloading = Button('Start/Stop autoreloading')

    def _start_stop_reloading_fired(self):
        if self.reloading_thread and self.reloading_thread.is_alive():
            self.reloading_thread.event.set()
        else:
            self.reloading_thread = ReloadingThread()
            #self.reloading_thread.time_to_sleep = self.time_to_sleep
            self.reloading_thread.controler = self
            self.reloading_thread.event = self.event
            self.reloading_thread.event.clear()
            self.reloading_thread.start()

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
        lines = []
        for lds in self.ldfiles.ldfiles:
            for ld in lds.ld_curves:
                lines.extend(ld.draw_line(ax))
        cursor = mplcursors.cursor(lines, highlight=True)
        cursor.connect("add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))
        if self.figure_settings.show_legend:
            ax.legend()
        if self.figure_settings.set_xlim:
            ax.set_xlim(self.figure_settings.x_limits)
        if self.figure_settings.set_ylim:
            ax.set_ylim(self.figure_settings.y_limits)
        self.figure.canvas.draw()
        print('redrawn all')

    def _reload_button_fired(self):
        for lds in self.ldfiles.ldfiles:
            lds.reload_button = True
        print('reloaded all')
        
    def _get_current_limits_fired(self):
        ax = self.figure.axes[0]
        self.figure_settings.x_limits = ax.get_xlim()
        self.figure_settings.y_limits = ax.get_ylim()

    view = View(Item('ldfiles', editor=tree_editor, show_label=False, id='treeeditor'),
                '_',
                Item('@figure_settings', show_label=False),
                Item('get_current_limits', show_label=False),
                HGroup(
                       Item('reload_button', show_label=False, springy=True),
                       Item('clear_button', show_label=False, springy=True),
                       Item('draw_button', show_label=False, springy=True)),
                HGroup(Item('time_to_sleep', enabled_when='reloading_thread.event.is_set()'),
                       Item('start_stop_reloading', show_label=False)),
                resizable=True,
                id='control_panel',
            )


class TC_Handler(Handler):

    def object_title_changed(self, info):
        if info.initialized:
            info.ui.title = info.object.title
            
    def close(self, info, is_OK):
        if ( info.object.panel.reloading_thread
            and info.object.panel.reloading_thread.is_alive() ):
            info.object.panel.reloading_thread.event.set()
            while info.object.panel.reloading_thread.is_alive():
                time.sleep(0.1)
        return True


class LDViewer(HasStrictTraits):

    title = Str('LD_Viewer')

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
    
    view = View(Item('title', label='window title'),
                HSplit('@panel',
                       Item('figure', show_label=False, editor=MPLFigureEditor(), dock='vertical', width=0.8),
                       show_labels=False, id='ld_hsplit', dock='tab'
                      ),
                title='LD Viewer',
                resizable=True,
                height=0.75, width=0.75,
                buttons=NoButtons,
                handler=TC_Handler(),
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
    start_time = time.time()
    logging.info(f'Loading data')
    if args.ld_files:
        for ld_idx, ld_file in enumerate(args.ld_files):
            ld_viewer.panel.ldfiles.add_ldfile = True
            ld_viewer.panel.ldfiles.ldfiles[ld_idx].ld_file = str(pathlib.Path(ld_file).absolute())
            ld_viewer.panel.ldfiles.ldfiles[ld_idx].name = pathlib.Path(ld_file).absolute().parts[-3]
            ld_viewer.panel.ldfiles.ldfiles[ld_idx].add_ldcurve = True
            if args.x:# in ld_viewer.labels:
                #ld_viewer.x_values = args.x
                ld_viewer.panel.ldfiles.ldfiles[ld_idx].ld_curves[0].x_values = args.x
            if args.y:# in ld_viewer.labels:
                #ld_viewer.y_values = args.y
                ld_viewer.panel.ldfiles.ldfiles[ld_idx].ld_curves[0].y_values = args.y
    logging.info(f'Data loaded {time.time() - start_time}')
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
