import os
import time
import pathlib
import argparse
import logging
import datetime
#os.environ["ETS_TOOLKIT"] = "qt"
#os.environ["QT_API"] = "pyqt6"
from traits.api import HasStrictTraits, Instance, List, Enum, Button, Str, \
                        Any, Bool, Float, Tuple
from traitsui.api import Item, View, HSplit, NoButtons, HGroup, Group, \
                         Handler, VGroup, ListEditor, TreeEditor, TreeNode
from traitsui.menu \
    import Menu, Action, Separator
from traitsui.qt.tree_editor \
    import CopyAction, CutAction, \
           PasteAction, DeleteAction, RenameAction
from traitsui.message import Message
from traitsui.file_dialog  import open_file, TextInfo, FileInfo, save_file
from mpl_qt_editor import MPLFigureEditor
from matplotlib.figure import Figure
import pandas as pd
import mplcursors
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


def error(message="", title="Error", parent=None):
    """Displays an error message to the user as a modal window."""
    msg = Message(message=message)
    msg.edit_traits(
        parent=parent,
        view=View(
            ["message~", "|<>"], 
            title=title, 
            buttons=["OK"], 
            kind="modal",
            resizable=True,
            scrollable=True
        ),
    )


class LDFile(HasStrictTraits):
    ldfile_num = 0

    def __init__(self, **traits):
        LDFile.ldfile_num += 1
        super().__init__(**traits)

    name = Str
    ld_file = Str('', changed=True, enter_set=True, auto_set=False)
    open_button = Button('Open...')
    reload_button = Button('Reload')
    data = Any(None)
    labels = List([], changed=True)

    add_ldcurve = Button('Add new LD curve')
    ld_curves = List
    figure = Instance(Figure)
    parent = Any()  # Reference to parent LDFiles object

    def _add_ldcurve_fired(self):
        if self.figure is None:
            logging.error("Figure not initialized in LDFile, cannot add LD curve")
            error("Figure not initialized. Please restart the application.", "Error")
            return
        self.ld_curves.append(LDCurve(ldfile=self, labels=self.labels))
    
    def duplicate(self):
        """Duplicate this LDFile with all its curves"""
        if self.parent is None:
            logging.error("Cannot duplicate: no parent reference")
            return
        
        new_file = LDFile(figure=self.figure, parent=self.parent)
        new_file.name = f"{self.name}_copy"
        new_file.ld_file = self.ld_file
        new_file.data = self.data
        new_file.labels = self.labels.copy() if self.labels else []
        
        # Duplicate all curves
        for curve in self.ld_curves:
            new_curve = LDCurve(ldfile=new_file, labels=new_file.labels)
            new_curve.name = curve.name
            new_curve.x_values = curve.x_values
            new_curve.y_values = curve.y_values
            new_curve.x_scale = curve.x_scale
            new_curve.y_scale = curve.y_scale
            new_file.ld_curves.append(new_curve)
        
        self.parent.ldfiles.append(new_file)
        logging.info(f'Duplicated {self.name} with {len(new_file.ld_curves)} curves')

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
        """Load data from LD file"""
        if not self.ld_file:
            logging.debug('No file specified yet')
            self.data = None
            return []
        
        if not pathlib.Path(self.ld_file).exists():
            logging.error(f'File does not exist: {self.ld_file}')
            error(f'File does not exist:\n{self.ld_file}', 'Error')
            self.data = None
            return []
        
        try:
            logging.debug(f'Loading: {self.ld_file}')
            self.data = pd.read_csv(self.ld_file, sep='\t', header=0, index_col=False)
            logging.info(f'Loaded {len(self.data)} rows from {pathlib.Path(self.ld_file).name}')
            return list(self.data.columns)
        except Exception as err:
            logging.error(f'Failed to load {self.ld_file}: {err}')
            error(f'The file cannot be loaded:\n{self.ld_file}\n\n{str(err)}', 'Error')
            self.data = None
            return []

    def _labels_changed(self):
        for curve in self.ld_curves:
            curve.labels = self.labels

    view = View(VGroup(HGroup(Item('open_button', show_label=False, id='ld_open'),
                       Item('reload_button', show_label=False),
                       Item('ld_file', id='ld_file'),
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

    draw_this_button = Button('Draw this one')
    copy_to_clipboard = Button('Copy to clipboard')
    save_to_txt = Button('Save to txt')

    @property
    def figure(self):
        """Get figure from parent LDFile"""
        return self.ldfile.figure if self.ldfile else None

    def _draw_this_button_fired(self):
        if not self._validate_figure():
            return
        ax = self.figure.axes[0]
        self._draw_this(ax)
        self.figure.canvas.draw()

    def _validate_figure(self):
        """Validate that figure is properly initialized"""
        if self.figure is None:
            logging.error("Figure not initialized")
            error("Figure not initialized. Cannot draw curve.", "Error")
            return False
        if not self.figure.axes:
            logging.error("Figure has no axes")
            error("Figure has no axes. Cannot draw curve.", "Error")
            return False
        return True

    def _draw_this(self, ax):
        """Draw this curve on the given axes"""
        self.draw_line(ax)

    def _validate_data(self):
        """Validate that data is available and contains required columns"""
        if self.ldfile.data is None:
            logging.warning(f"No data available for {self.ldfile.name}-{self.name}")
            return False
        if not isinstance(self.ldfile.data, pd.DataFrame) or len(self.ldfile.data) == 0:
            logging.warning(f"Invalid or empty data for {self.ldfile.name}-{self.name}")
            return False
        if self.x_values not in self.ldfile.data.columns or self.y_values not in self.ldfile.data.columns:
            logging.error(f"Column {self.x_values} or {self.y_values} not found in data")
            return False
        return True

    def draw_line(self, ax):
        """Draw line on the given axes"""
        if not self._validate_data():
            return []
        
        x_to_draw = self.ldfile.data[self.x_values] * self.x_scale
        y_to_draw = self.ldfile.data[self.y_values] * self.y_scale
        
        return ax.plot(x_to_draw, y_to_draw,
                       label=f'{self.ldfile.name}-{self.name}', 
                       marker='o')

    def _copy_to_clipboard_fired(self):
        if not self._validate_data():
            error("No data available to copy.", "Error")
            return
        
        try:
            x_data = self.ldfile.data[self.x_values] * self.x_scale
            y_data = self.ldfile.data[self.y_values] * self.y_scale
            df = pd.DataFrame({self.x_values: x_data, self.y_values: y_data})
            df.to_clipboard(excel=True, index=False)
            logging.info('Data copied to clipboard')
        except Exception as err:
            logging.error(f"Clipboard copy failed: {err}")
            error(f"Failed to copy to clipboard: {str(err)}", "Error")

    def _save_to_txt_fired(self):
        if not self._validate_data():
            error("No data available to save.", "Error")
            return
        
        file_name = save_file(extensions=FileInfo(), id='savefile')
        if not file_name:
            logging.info('Save cancelled - no filename provided.')
            return
        
        try:
            x_data = self.ldfile.data[self.x_values] * self.x_scale
            y_data = self.ldfile.data[self.y_values] * self.y_scale
            df = pd.DataFrame({self.x_values: x_data, self.y_values: y_data})
            df.to_csv(file_name, sep='\t', index=False)
            logging.info(f'Data saved to {file_name}')
        except Exception as err:
            logging.error(f"File save failed: {err}")
            error(f"Failed to save file: {str(err)}", "Error")

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
    ldfiles = List(LDFile)

    def _add_ldfile_fired(self):
        if self.figure is None:
            logging.error("Figure not initialized, cannot add LD file")
            error("Figure not initialized. Please restart the application.", "Error")
            return
        new_file = LDFile(figure=self.figure, parent=self)
        self.ldfiles.append(new_file)
    
    def _ldfiles_items_changed(self, event):
        """Ensure parent is set when items are added to ldfiles list"""
        for item in event.added:
            if item.parent is None:
                item.parent = self
            if item.figure is None:
                item.figure = self.figure
    
    def duplicate_ldfile(self, ldfile):
        """Duplicate an LDFile with all its curves"""
        if not isinstance(ldfile, LDFile):
            return
        
        new_file = LDFile(figure=self.figure)
        new_file.name = f"{ldfile.name}_copy"
        new_file.ld_file = ldfile.ld_file
        new_file.data = ldfile.data
        new_file.labels = ldfile.labels.copy() if ldfile.labels else []
        
        # Duplicate all curves
        for curve in ldfile.ld_curves:
            new_curve = LDCurve(ldfile=new_file, labels=new_file.labels)
            new_curve.name = curve.name
            new_curve.x_values = curve.x_values
            new_curve.y_values = curve.y_values
            new_curve.x_scale = curve.x_scale
            new_curve.y_scale = curve.y_scale
            new_file.ld_curves.append(new_curve)
        
        self.ldfiles.append(new_file)
        logging.info(f'Duplicated {ldfile.name} with {len(new_file.ld_curves)} curves')

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
    x_limits = Tuple(float('nan'), float('nan'))
    set_ylim = Bool()
    y_limits = Tuple(float('nan'), float('nan'))
    view = View(
        Item('ignore_last_step', show_label=True),
        Item('show_legend', show_label=True),
        HGroup(Item('set_xlim', show_label=False),
               Item('x_limits', show_label=True)),
        HGroup(Item('set_ylim', show_label=False),
               Item('y_limits', show_label=True)),
        resizable=True
    )


class AxisSettings(HasStrictTraits):
    title = Str('')
    xlabel = Str('')
    ylabel = Str('')
    
    view = View(
        Item('title', show_label=True),
        Item('xlabel', show_label=True),
        Item('ylabel', show_label=True),
        resizable=True
    )
        

class ReloadingThread(Thread):
    """Background thread for auto-reloading data at intervals"""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.daemon = True  # Thread dies when main program exits

    def run(self):
        while not self.event.is_set():
            try:
                self.controler._reload_button_fired()
                self.controler._clear_button_fired()
                self.controler._draw_button_fired()
            except Exception as e:
                logging.error(f"Error in reloading thread: {e}")
            self.event.wait(timeout=self.controler.time_to_sleep)



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
                
# Custom action for adding new LD file with proper figure initialization
add_ldfile_action = Action(
    name='New LD File',
    action='object._add_ldfile_fired'
)

# Custom action for duplicating an LD file
duplicate_ldfile_action = Action(
    name='Duplicate',
    action='object.duplicate'
)

# Tree editor
tree_editor = TreeEditor(
    nodes = [
        TreeNode( node_for  = [ LDFiles ],
                  auto_open = True,
                  children = 'ldfiles',
                  label     = '=LD Files',
                  add       = [ LDFile ],
                  move      = [ LDFile ],
                  menu=Menu( add_ldfile_action,
                             Separator(),
                             CopyAction,
                             CutAction,
                             PasteAction,
                             Separator(),
                             DeleteAction),
                  view =  View()),
        TreeNode( node_for  = [ LDFile ],
                  auto_open = True,
                  label     = 'name',
                  copy       = True,
                  menu=Menu( add_ldfile_action,
                             Separator(),
                             CopyAction,
                             duplicate_ldfile_action,
                             CutAction,
                             PasteAction,
                             Separator(),
                             DeleteAction,
                             Separator(),
                             RenameAction),
                  view =  ldfile_view)
    ],
    hide_root=False,
    orientation="vertical",
)
                  

class ControlPanel(HasStrictTraits):
    ldfiles = Instance(LDFiles)
    figure_settings = Instance(FigureSettings, ())
    axis_settings = Instance(AxisSettings, ())
    figure = Instance(Figure)
    cursor = Any()  # mplcursors cursor instance

    add_ldfile_button = Button('Add new LD file')
    draw_button = Button('Draw all')
    reload_button = Button('Reload all')
    clear_button = Button('Clear')
    update_axis_labels = Button('Update default labels')
    
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

    def _ldfiles_default(self):
        logging.debug(f'Initializing LDFiles with figure: {self.figure}')
        ldfiles = LDFiles(figure=self.figure)
        # Create one empty LDFile on startup
        ldfiles.ldfiles.append(LDFile(figure=self.figure, parent=ldfiles))
        return ldfiles

    def _add_ldfile_button_fired(self):
        self.ldfiles.add_ldfile = True
    
    def _update_axis_labels_fired(self):
        """Generate default axis labels from all curves"""
        x_values_set = set()
        y_values_set = set()
        
        for ldfile in self.ldfiles.ldfiles:
            for curve in ldfile.ld_curves:
                if curve.x_values:
                    x_values_set.add(curve.x_values)
                if curve.y_values:
                    y_values_set.add(curve.y_values)
        
        if x_values_set:
            self.axis_settings.xlabel = ', '.join(sorted(x_values_set))
        if y_values_set:
            self.axis_settings.ylabel = ', '.join(sorted(y_values_set))
        
        logging.info(f'Updated axis labels: x="{self.axis_settings.xlabel}", y="{self.axis_settings.ylabel}"')

    def _clear_button_fired(self):
        # Remove old cursor if exists
        if self.cursor is not None:
            self.cursor.remove()
            self.cursor = None
        ax = self.figure.axes[0]
        ax.cla()
        self.figure.canvas.draw()

    def _draw_button_fired(self):
        if not self.figure or not self.figure.axes:
            logging.error('Figure not properly initialized')
            error('Figure not properly initialized', 'Error')
            return
        
        # Remove old cursor to clear any active annotations
        if self.cursor is not None:
            self.cursor.remove()
            self.cursor = None
        
        ax = self.figure.axes[0]
        lines = []
        
        for ldfile in self.ldfiles.ldfiles:
            for curve in ldfile.ld_curves:
                line_result = curve.draw_line(ax)
                if line_result:
                    lines.extend(line_result)
        
        # Create cursor for the new lines
        if lines:
            self.cursor = mplcursors.cursor(lines, highlight=True)
            self.cursor.connect("add", lambda sel: sel.annotation.set_text(sel.artist.get_label()))
        
        # Apply axis settings
        if self.axis_settings.title:
            ax.set_title(self.axis_settings.title)
        if self.axis_settings.xlabel:
            ax.set_xlabel(self.axis_settings.xlabel)
        if self.axis_settings.ylabel:
            ax.set_ylabel(self.axis_settings.ylabel)
        
        if self.figure_settings.show_legend:
            ax.legend()
        if self.figure_settings.set_xlim:
            ax.set_xlim(self.figure_settings.x_limits)
        if self.figure_settings.set_ylim:
            ax.set_ylim(self.figure_settings.y_limits)
        
        self.figure.canvas.draw()
        logging.info(f'Drew {len(lines)} curve(s)')

    def _reload_button_fired(self):
        for lds in self.ldfiles.ldfiles:
            lds.reload_button = True
        logging.info('Reloaded all files')
        
    def _get_current_limits_fired(self):
        ax = self.figure.axes[0]
        self.figure_settings.x_limits = ax.get_xlim()
        self.figure_settings.y_limits = ax.get_ylim()

    view = View(VGroup(
                    Group(
                        Group(
                            Item('add_ldfile_button', show_label=False),
                            Item('ldfiles', editor=tree_editor, show_label=False, id='treeeditor'), label='LD Files',
                            id='ld_files_tab'),
                        Group(
                            Item('@figure_settings', show_label=False),
                            Item('get_current_limits', show_label=False),
                            '_',
                            Item('@axis_settings', show_label=False),
                            Item('update_axis_labels', show_label=False), 
                            label='Figure Settings',
                            id='figure_settings_tab'),
                        layout='tabbed',
                        id='ld_control_tabs',
                    ),
                    HGroup(
                        Item('reload_button', show_label=False, springy=True),
                        Item('clear_button', show_label=False, springy=True),
                        Item('draw_button', show_label=False, springy=True)),
                    HGroup(Item('time_to_sleep', enabled_when='reloading_thread.event.is_set()'),
                        Item('start_stop_reloading', show_label=False)),),
                resizable=True,
                id='control_panel',
            )


class TC_Handler(Handler):

    def object_title_changed(self, info):
        if info.initialized:
            info.ui.title = info.object.title
            
    def close(self, info, is_OK):
        """Clean up before closing"""
        if (info.object.panel.reloading_thread
            and info.object.panel.reloading_thread.is_alive()):
            info.object.panel.reloading_thread.event.set()
            # Wait briefly for thread to finish (it's daemon so will die anyway)
            info.object.panel.reloading_thread.join(timeout=0.5)
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
    ld_viewer = LDViewer()
    start_time = time.time()
    logging.info('Loading data')
    
    if args.ld_files:
        for ld_idx, ld_file in enumerate(args.ld_files):
            ld_viewer.panel.ldfiles.add_ldfile = True
            ld_viewer.panel.ldfiles.ldfiles[ld_idx].ld_file = str(pathlib.Path(ld_file).absolute())
            ld_viewer.panel.ldfiles.ldfiles[ld_idx].name = pathlib.Path(ld_file).absolute().parts[-3]
            ld_viewer.panel.ldfiles.ldfiles[ld_idx].add_ldcurve = True
            
            if args.x:
                ld_viewer.panel.ldfiles.ldfiles[ld_idx].ld_curves[0].x_values = args.x
            if args.y:
                ld_viewer.panel.ldfiles.ldfiles[ld_idx].ld_curves[0].y_values = args.y
    
    logging.info(f'Data loaded in {time.time() - start_time:.2f}s')
    ld_viewer.configure_traits()
