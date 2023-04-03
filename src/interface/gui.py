import base64
from dataclasses import dataclass
from io import BytesIO
import os

from copy import deepcopy
from pathlib import Path
from threading import Thread
from typing import Any

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
    self.__window: sg.Window = None              # main frame window
    self.__results_folder = results_folder       # where to find the system tree
    self.__font = ('Helvetica', 11)              # font for the main window
    self.__monospace_font = ('Courier', 7)       # font for the log window
    self.__file_tree: list[list[sg.Tree]] = None # file tree layout
    self.__tree_data: sg.TreeData = None         # file tree data

    self.__selected_organisms: dict[str, int] = {} # selected organisms
    self.__selected_region: str = 'CDS'            # selected region
    self.__tree = tree                             # tree stored to recieve `Value` data

    # this one bellow is stored here
    # because the reset thread recieve `self` as argument

    self.__tree_component = None # tree component (to update it)
    self.__data_component: dict[str, dict[str, Any]] = {}

    self.__known_files: set[str] = set() # known files

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

    path = self.__results_folder if path is None else path
    self.__tree_data = sg.TreeData()

    self.__tree_data.insert('', 0, os.path.basename(self.__results_folder), [], icon=folder_icon)

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
                sg.Button('Test', font=self.__font, size=(10, 2)),
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
        data=sg.TreeData(),      # empty tree, will be filled later
        headings=[],             # no headings
        auto_size_columns=False,
        num_rows=25,
        col0_width=40,
        key='-TREE-',
        show_expanded=True,
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

    capture.enable_proxy() # enable capture stdout and stderr for all children threads

    def new_key() -> int: # get a new key for the tree
      key = 1
      while key in self.__data_component:
        key += 1
      return key

    def dry_reload(gui: 'GenomeGUI') -> None: # reload the tree without reseting the data
      gui.__known_files.clear()
      gui.__data_component = deepcopy(base_data)
      gui.__tree_component.update(values=gui.__build_tree())

    def reset(gui: 'GenomeGUI') -> None: # reset the system tree
      gui.__tree.build(True, True)
      load_tree(gui)
      gui.__known_files.clear()
      gui.__data_component = deepcopy(base_data)
      gui.__tree_component.update(values=gui.__tree_data)

    def get_data(gui: 'GenomeGUI', val: Value, parent_key: int) -> None:
      create_data_from_NC(val.name, val.path, val.nc, gui.__selected_region)
      files = os.listdir(val.path)
      for f in files:
        fullname = os.path.join(val.path, f)
        if fullname not in gui.__known_files:
          gui.__known_files.add(f)
          key = new_key()
          gui.__data_component[key] = {
            'kind': FILE,
            'path': val.path,
            'file': f,
            'children': None,
          }
          gui.__tree_data.insert(parent_key, key, f, [], icon=file_icon)
          gui.__known_files.add(f)
      gui.__tree_component.update(values=gui.__tree_data)
      iid = gui.__tree_component.KeyToID[parent_key]
      gui.__tree_component.Widget.see(iid)

    self.__tree_component = self.__window['-TREE-'] # get the tree component
    if isinstance(self.__tree_component, sg.ErrorElement):
      panic('Tree not found')

    self.__tree_component.Widget.configure(show='tree')
    self.__tree_component.bind('<Double-1>', "DOUBLE-CLICK-") # expand/explore on double click
    self.__tree_component.bind('<Button-3>', "RIGHT-CLICK-")  # select on right click

    def load_tree(gui: 'GenomeGUI') -> None:
      gui.__tree_component.update(values=gui.__build_tree())
      info('Tree has been loaded')
      return

    Thread(target=load_tree, args=(self,)).start() # update the tree in a new thread
                                                   # load_tree(self)

    DIR, FILE = True, False
    self.__data_component = {
      0: {
        'kind': DIR,
        'path': '',
        'file': os.path.basename(self.__results_folder),
        'children': None
      },
    }
    base_data = deepcopy(self.__data_component) # base data to reset the tree

    running = True # main loop
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
        keys: list[int] = []
        for organism in self.__selected_organisms.keys():
          try:
            values.append(self.__tree.get_info(organism))
            keys.append(self.__selected_organisms[organism])
          except KeyError:
            error(f'No data or key for {organism} : '
                  f'this could be because {organism} is not a valid organism')
            if len(values) != len(keys):
              # here normally we only have one item more in values
              # because we failed to get the key
              values.pop()
            continue
        for i, val in enumerate(values):
          Thread(target=get_data, args=(self, val, keys[i])).start()
        for organism, val in self.__selected_organisms.items():
          full_path = self.__data_component[val]['path']
          self.__tree_component.update(key=val, icon=folder_icon if os.path.isdir(full_path) else file_icon)
        self.__selected_organisms.clear()

      if event == 'Dry Reload':
        Thread(target=dry_reload, args=(self,)).start()
      if event == 'Reset':
        r = sg.popup_ok_cancel('Are you sure you want to reset the tree?\n'
                               'This will delete all data in "Results" folder!')
        if r == 'OK':
          Thread(target=reset, args=(self,)).start()

      if event in {'-TREE-RIGHT-CLICK-', 'Toggle Selection'}:
        for selected in values['-TREE-']:
          name = self.__data_component[selected]['file']
          full_path = self.__data_component[selected]['path']
          full_path = os.path.join(full_path, name)
          if name in self.__selected_organisms.keys():
            self.__selected_organisms.pop(name)
            self.__tree_component.update(key=selected,
                                         icon=folder_icon if os.path.isdir(full_path) else file_icon)
          else:
            self.__selected_organisms[name] = selected
            self.__tree_component.update(
              key=selected, icon=folder_checked_icon if os.path.isdir(full_path) else file_checked_icon)
          info(f'Selected organisms: {self.__selected_organisms.keys()}')

      if event == '-TREE-DOUBLE-CLICK-':
        try:
          parent_key = values['-TREE-'][0]
          node = self.__data_component[parent_key]
        except IndexError:
          continue

        if node['kind'] == DIR:
          parent_path = Path(node['path']).joinpath(node['file'])
          files: list
          try:
            files = sorted(list(parent_path.iterdir()), key=lambda file: file.is_file())
          except PermissionError:
            error(f'Permission denied: {parent_path}')
            continue

          if node['children'] is not None:
            # here we only search for new files
            for item in files:
              if item.is_file() and str(item) not in self.__known_files:
                self.__known_files.add(str(item))
                key = new_key()
                kind, path, file = item.is_dir(), str(item.parent), item.name
                self.__tree_data.insert(parent_key,
                                        key,
                                        str(file), [],
                                        icon=folder_icon if kind == DIR else file_icon)
                node['children'].append(key)
                self.__data_component[key] = {'kind': kind, 'path': path, 'file': file, 'children': None}
          else:
            # here we search for new files and new folders that were never expanded
            node['children'] = []
            for item in files:
              kind, path, file = item.is_dir(), str(item.parent), item.name
              if kind == FILE and str(item.name) in self.__known_files:
                continue
              key = new_key()
              self.__known_files.add(str(item))
              self.__tree_data.insert(parent_key,
                                      key,
                                      str(file), [],
                                      icon=folder_icon if kind == DIR else file_icon)
              node['children'].append(key)
              self.__data_component[key] = {'kind': kind, 'path': path, 'file': file, 'children': None}

          debug(f'Loaded {len(files)} items from {parent_path}')
          self.__tree_component.update(values=self.__tree_data)
          iid = self.__tree_component.KeyToID[parent_key]
          self.__tree_component.Widget.see(iid)
        
        else:
          # here we open the file
          path = Path(node['path']).joinpath(node['file'])
          try:
            info(f'Opening {path}')
            os.startfile(path, 'open')
          except PermissionError:
            error(f'Permission denied: {path}')

    capture.disable_proxy()
    self.__window.close()
