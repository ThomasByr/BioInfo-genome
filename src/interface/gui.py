import base64
from io import BytesIO
import os

from threading import Thread
from multiprocessing import Pool

from PIL import Image, ImageDraw
import PySimpleGUI as sg

from ..core import Tree, Value, create_data_from_NC
from ..helper import debug, info, error, panic, capture

__all__ = ['GenomeGUI']

folder_icon = b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAACXBIWXMAAAsSAAALEgHS3X78AAABnUlEQVQ4y8WSv2rUQRSFv7vZgJFFsQg2EkWb4AvEJ8hqKVilSmFn3iNvIAp21oIW9haihBRKiqwElMVsIJjNrprsOr/5dyzml3UhEQIWHhjmcpn7zblw4B9lJ8Xag9mlmQb3AJzX3tOX8Tngzg349q7t5xcfzpKGhOFHnjx+9qLTzW8wsmFTL2Gzk7Y2O/k9kCbtwUZbV+Zvo8Md3PALrjoiqsKSR9ljpAJpwOsNtlfXfRvoNU8Arr/NsVo0ry5z4dZN5hoGqEzYDChBOoKwS/vSq0XW3y5NAI/uN1cvLqzQur4MCpBGEEd1PQDfQ74HYR+LfeQOAOYAmgAmbly+dgfid5CHPIKqC74L8RDyGPIYy7+QQjFWa7ICsQ8SpB/IfcJSDVMAJUwJkYDMNOEPIBxA/gnuMyYPijXAI3lMse7FGnIKsIuqrxgRSeXOoYZUCI8pIKW/OHA7kD2YYcpAKgM5ABXk4qSsdJaDOMCsgTIYAlL5TQFTyUIZDmev0N/bnwqnylEBQS45UKnHx/lUlFvA3fo+jwR8ALb47/oNma38cuqiJ9AAAAAASUVORK5CYII='
file_icon = b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAACXBIWXMAAAsSAAALEgHS3X78AAABU0lEQVQ4y52TzStEURiHn/ecc6XG54JSdlMkNhYWsiILS0lsJaUsLW2Mv8CfIDtr2VtbY4GUEvmIZnKbZsY977Uwt2HcyW1+dTZvt6fn9557BGB+aaNQKBR2ifkbgWR+cX13ubO1svz++niVTA1ArDHDg91UahHFsMxbKWycYsjze4muTsP64vT43v7hSf/A0FgdjQPQWAmco68nB+T+SFSqNUQgcIbN1bn8Z3RwvL22MAvcu8TACFgrpMVZ4aUYcn77BMDkxGgemAGOHIBXxRjBWZMKoCPA2h6qEUSRR2MF6GxUUMUaIUgBCNTnAcm3H2G5YQfgvccYIXAtDH7FoKq/AaqKlbrBj2trFVXfBPAea4SOIIsBeN9kkCwxsNkAqRWy7+B7Z00G3xVc2wZeMSI4S7sVYkSk5Z/4PyBWROqvox3A28PN2cjUwinQC9QyckKALxj4kv2auK0xAAAAAElFTkSuQmCC'

folder_img = Image.open(BytesIO(base64.b64decode(folder_icon)))
file_img = Image.open(BytesIO(base64.b64decode(file_icon)))
tick_lines = [(4, 8), (7, 11), (12, 4)]

folder_draw = ImageDraw.Draw(folder_img)
folder_draw.line(tick_lines, fill=(0, 155, 0), width=2)
with BytesIO() as output:
  folder_img.save(output, format='PNG')
  folder_checked_icon = output.getvalue()

file_draw = ImageDraw.Draw(file_img)
file_draw.line(tick_lines, fill=(0, 155, 0), width=2)
with BytesIO() as output:
  file_img.save(output, format='PNG')
  file_checked_icon = output.getvalue()


class GenomeGUI:

  def __init__(self, tree: Tree, results_folder: str = 'Results'):
    self.__window: sg.Window = None  # main frame window
    self.__results_folder = results_folder  # where to find the system tree
    self.__font = ('Helvetica', 11)  # font for the main window
    self.__monospace_font = ('Courier', 7)  # font for the log window
    self.__file_tree: list[list[sg.Tree]] = None  # file tree layout
    self.__tree_data: sg.TreeData = None  # file tree data

    self.__selected_organisms: set[str] = set()  # selected organisms
    self.__selected_region: str = 'CDS'  # selected region
    self.__tree = tree  # tree stored to recieve `Value` data

    # this one bellow is stored here
    # because the reset thread recieve `self` as argument

    self.__tree_component = None

    # some layouts

    self.__info = [[sg.Text('Region Selection', font=self.__font)],
                   [
                     sg.Combo(
                       [
                         'CDS',
                         'centromere',
                         'intron',
                         'mobile_element',
                         'ncRNA',
                         'rRNA',
                         'telomere',
                         'tRNA',
                         '3\'UTR',
                         '5\'UTR',
                       ],
                       enable_events=True,
                       expand_x=True,
                       default_value='CDS',
                       key='-SELECT-',
                     )
                   ]]

    self.__log = [[sg.Text('Logs', font=self.__font)],
                  [
                    sg.Multiline(
                      key='-LOG-',
                      size=(60, 25),
                      expand_x=True,
                      echo_stdout_stderr=True,
                      reroute_stdout=True,
                      reroute_stderr=True,
                      autoscroll=True,
                      write_only=True,
                      font=self.__monospace_font,
                    )
                  ]]
    self.__right_side = [
      [sg.Column(self.__info, vertical_alignment='top')],
      [sg.HSeparator()],
      [sg.Column(self.__log, vertical_alignment='top')],
    ]

  def __build_tree(self, path: str = None) -> sg.TreeData:

    def add_folder(parent: str, folder: str):
      files = os.listdir(folder)
      for f in files:
        fullname = os.path.join(folder, f)
        if os.path.isdir(fullname):
          self.__tree_data.insert(parent, fullname, f, values=[], icon=folder_icon)
          add_folder(fullname, fullname)
        # else:
        #   self.__tree_data.insert(parent, fullname, f, values=[], icon=file_icon)

    path = self.__results_folder if path is None else path
    self.__tree_data = sg.TreeData()

    add_folder('', path)

    return self.__tree_data

  def __build_layout(self):
    layout = [[sg.Text('GENOME', font=self.__font)], [sg.HSeparator()],
              [
                sg.Column(self.__file_tree, justification='center'),
                sg.VSeparator(),
                sg.Column(self.__right_side)
              ], [sg.HSeparator()],
              [
                sg.Button('Run', font=self.__font, size=(10, 2)),
                sg.Button('Toggle Selection', font=self.__font, size=(10, 2)),
                sg.VSeparator(),
                sg.Button('Dry Reload', font=self.__font, size=(10, 2)),
                sg.Button('Reset', font=self.__font, size=(10, 2)),
              ], [sg.Text('Developed with ❤', font=self.__font, justification='right', expand_x=True)],
              [
                sg.Text('by: @ThomasByr, @m7415, @JBrandstaedt and @Bas6700',
                        font=self.__font,
                        justification='right',
                        expand_x=True)
              ]]
    return layout

  def __build_file_tree(self):
    self.__file_tree = [[
      sg.Tree(
        data=sg.TreeData(),  # empty tree, will be filled later
        headings=[],
        auto_size_columns=False,
        num_rows=25,
        col0_width=40,
        key='-TREE-',
        show_expanded=False,
        enable_events=True,
      )
    ]]

  def __build_window(self):
    bioinformatics: bytes
    with open(os.path.join('assets', 'bioinformatics.png'), 'rb') as f:
      bioinformatics = f.read()
    self.__window = sg.Window(
      'Genome',
      self.__build_layout(),
      size=(800, 600),
      auto_size_text=False,
      auto_size_buttons=False,
      icon=bioinformatics,
      finalize=True,
    )

  def run(self):
    """
    run the GUI\n
    captures stdout and stderr and quits on window close
    """
    sg.theme('DarkTeal9')
    self.__build_file_tree()
    self.__build_window()

    capture.enable_proxy()  # enable capture stdout and stderr for all children threads

    def reset(gui: 'GenomeGUI') -> None:  # reset the system tree
      gui.__tree.build(True, True)
      gui.__build_file_tree()
      gui.__tree_component.update(values=gui.__tree_data)

    def get_data(gui: 'GenomeGUI', val: Value) -> None:
      create_data_from_NC(val.name, val.path, val.nc, gui.__selected_region)

    self.__tree_component = self.__window['-TREE-']  # get the tree component
    if isinstance(self.__tree_component, sg.ErrorElement):
      panic('Tree not found')

    self.__tree_component.Widget.configure(show='tree')
    self.__tree_component.bind('<Button-3>', "RIGHT-CLICK-")  # select on right click

    def load_tree(gui: 'GenomeGUI') -> None:
      gui.__tree_component.update(values=gui.__build_tree())
      info('Tree has been loaded')
      return

    Thread(target=load_tree, args=(self,)).start()  # update the tree in a new thread

    running = True  # main loop
    key = None
    while running:
      event, values = self.__window.read()
      if event in (sg.WIN_CLOSED, 'Exit'):
        running = False

      if event == '-SELECT-':
        self.__selected_region = values['-SELECT-']
        info(f'Selected region: {self.__selected_region}')

      if event == 'Run':
        if len(self.__selected_organisms) == 0:
          error('No organism selected')
          continue
        if self.__selected_region is None:
          error('No region selected')
          continue
        values: list[Value] = []
        for organism in self.__selected_organisms:
          try:
            values.append(self.__tree.get_info(organism))
          except KeyError:
            error(f'No data for {organism} : '
                  f'this could be because {organism} is not a valid organism')
            continue
        for val in values:
          Thread(target=get_data, args=(self, val)).start()

      if event == 'Dry Reload':
        Thread(target=self.__tree.build, args=[False, True]).start()
      if event == 'Reset':
        r = sg.popup_ok_cancel('Are you sure you want to reset the tree?\n'
                               'This will delete all data in "Results" folder!')
        if r == 'OK':
          Thread(target=reset, args=[self]).start()

      if event in {'-TREE-RIGHT-CLICK-', 'Toggle Selection'}:
        for item in self.__tree_component.Widget.selection():
          key = self.__tree_component.IdToKey[item]
          name = key.split(os.sep)[-1]
          if name in self.__selected_organisms:
            self.__selected_organisms.remove(name)
            self.__tree_component.update(key=key, icon=folder_icon if os.path.isdir(key) else file_icon)
          else:
            self.__selected_organisms.add(name)
            self.__tree_component.update(
              key=key, icon=folder_checked_icon if os.path.isdir(key) else file_checked_icon)
          info(f'Selected organisms: {self.__selected_organisms}')

    capture.disable_proxy()
    self.__window.close()
