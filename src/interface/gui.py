import os

from copy import deepcopy
from threading import Thread

from pathlib import Path
from typing import Any
import PySimpleGUI as sg

from ..core import Tree, Value, create_data_from_NC
from ..helper import debug, info, error, panic, capture

__all__ = ['GenomeGUI']

folder_icon = b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAACXBIWXMAAAsSAAALEgHS3X78AAABnUlEQVQ4y8WSv2rUQRSFv7vZgJFFsQg2EkWb4AvEJ8hqKVilSmFn3iNvIAp21oIW9haihBRKiqwElMVsIJjNrprsOr/5dyzml3UhEQIWHhjmcpn7zblw4B9lJ8Xag9mlmQb3AJzX3tOX8Tngzg349q7t5xcfzpKGhOFHnjx+9qLTzW8wsmFTL2Gzk7Y2O/k9kCbtwUZbV+Zvo8Md3PALrjoiqsKSR9ljpAJpwOsNtlfXfRvoNU8Arr/NsVo0ry5z4dZN5hoGqEzYDChBOoKwS/vSq0XW3y5NAI/uN1cvLqzQur4MCpBGEEd1PQDfQ74HYR+LfeQOAOYAmgAmbly+dgfid5CHPIKqC74L8RDyGPIYy7+QQjFWa7ICsQ8SpB/IfcJSDVMAJUwJkYDMNOEPIBxA/gnuMyYPijXAI3lMse7FGnIKsIuqrxgRSeXOoYZUCI8pIKW/OHA7kD2YYcpAKgM5ABXk4qSsdJaDOMCsgTIYAlL5TQFTyUIZDmev0N/bnwqnylEBQS45UKnHx/lUlFvA3fo+jwR8ALb47/oNma38cuqiJ9AAAAAASUVORK5CYII='
file_icon = b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAACXBIWXMAAAsSAAALEgHS3X78AAABU0lEQVQ4y52TzStEURiHn/ecc6XG54JSdlMkNhYWsiILS0lsJaUsLW2Mv8CfIDtr2VtbY4GUEvmIZnKbZsY977Uwt2HcyW1+dTZvt6fn9557BGB+aaNQKBR2ifkbgWR+cX13ubO1svz++niVTA1ArDHDg91UahHFsMxbKWycYsjze4muTsP64vT43v7hSf/A0FgdjQPQWAmco68nB+T+SFSqNUQgcIbN1bn8Z3RwvL22MAvcu8TACFgrpMVZ4aUYcn77BMDkxGgemAGOHIBXxRjBWZMKoCPA2h6qEUSRR2MF6GxUUMUaIUgBCNTnAcm3H2G5YQfgvccYIXAtDH7FoKq/AaqKlbrBj2trFVXfBPAea4SOIIsBeN9kkCwxsNkAqRWy7+B7Z00G3xVc2wZeMSI4S7sVYkSk5Z/4PyBWROqvox3A28PN2cjUwinQC9QyckKALxj4kv2auK0xAAAAAElFTkSuQmCC'


class GenomeGUI:

  def __init__(self, tree: Tree, results_folder: str = 'Results'):
    self.__window: sg.Window = None  # main frame window
    self.__results_folder = results_folder  # where to find the system tree
    self.__font = ('Helvetica', 11)  # font for the main window
    self.__monospace_font = ('Courier', 7)  # font for the log window
    self.__file_tree: list[list[sg.Tree]] = None  # file tree layout
    self.__tree_data: sg.TreeData = None  # file tree data
    self.__selected_organism: str = None  # selected organism
    self.__selected_region: str = 'CDS'  # selected region
    self.__tree = tree  # tree stored to recieve `Value` data

    # these two bellow are stored here
    # because the reset thread recieve `self` as argument

    self.__data_component: dict[str, dict[str, Any]] = {}
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
                  sg.Button('Confirm Selection', font=self.__font, size=(10, 2)),
                  sg.VSeparator(),
                  sg.Button('Dry Reload', font=self.__font, size=(10, 2)),
                  sg.Button('Reset', font=self.__font, size=(10, 2)),
              ], [sg.Text('Developed with ❤', font=self.__font, justification='right', expand_x=True)],
              [
                  sg.Text('by: @ThomasByr, @m7415, @JBrandstaedt and @Bas6700',
                          font=self.__font,
                          justification='right', expand_x=True)
              ]]
    return layout

  def __build_file_tree(self):
    self.__file_tree = [[
        sg.Tree(
            data=self.__build_tree(),
            headings=[],
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

    capture.enable_proxy()  # enable capture stdout and stderr for all children threads

    def new_key() -> int:  # get a new key for the tree
      key = 1
      while key in self.__data_component:
        key += 1
      return key

    def reset(gui: 'GenomeGUI') -> None:  # reset the system tree
      gui.__tree.build(True, True)
      gui.__build_file_tree()
      gui.__data_component = deepcopy(base_data)
      gui.__tree_component.update(values=gui.__tree_data)

    self.__tree_component = self.__window['-TREE-']  # get the tree component
    if isinstance(self.__tree_component, sg.ErrorElement):
      panic('Tree not found')

    self.__tree_component.Widget.configure(show='tree')
    self.__tree_component.bind('<Double-1>', "DOUBLE-CLICK-")  # expand/explore on double click

    DIR, FILE = True, False
    self.__data_component = {
        0: {
            'kind': DIR,
            'path': '',
            'file': os.path.basename(self.__results_folder),
            'children': None
        },
    }
    displayed_files: set[str] = set()  # do NOT expand a file twice
    base_data = deepcopy(self.__data_component)  # base data to reset the tree

    running = True  # main loop
    while running:
      event, values = self.__window.read()
      if event in (sg.WIN_CLOSED, 'Exit'):
        running = False

      if event == 'Confirm Selection':
        try:
          selected = values['-TREE-'][0]
          self.__selected_organism = self.__data_component[selected]['file']
          info(f'Selected organism: {self.__selected_organism}')
        except IndexError:
          error('Nothing selected in the tree')

      if event == '-SELECT-':
        self.__selected_region = values['-SELECT-']
        info(f'Selected region: {self.__selected_region}')

      if event == 'Run':
        if self.__selected_organism is None:
          error('No organism selected')
          continue
        if self.__selected_region is None:
          error('No region selected')
          continue
        val: Value
        try:
          val = self.__tree.get_info(self.__selected_organism)
        except KeyError:
          error(f'No data for {self.__selected_organism} : '
                f'this could be because {self.__selected_organism} is not a valid organism')
          continue
        Thread(target=create_data_from_NC,
               args=(self.__selected_organism, val.path, val.nc, self.__selected_region)).start()

      if event == 'Dry Reload':
        Thread(target=self.__tree.build, args=[False, True]).start()
      if event == 'Reset':
        r = sg.popup_ok_cancel('Are you sure you want to reset the tree?\n'
                               'This will delete all data in "Results" folder!')
        if r == 'OK':
          displayed_files.clear()
          Thread(target=reset, args=[self]).start()

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

          if node['children'] is not None:  # here we only search for new files
            for item in files:
              if item.is_file() and str(item) not in displayed_files:
                displayed_files.add(str(item))
                key = new_key()
                kind, path, file = item.is_dir(), str(item.parent), item.name
                self.__tree_data.insert(parent_key,
                                        key,
                                        str(file), [],
                                        icon=folder_icon if kind == DIR else file_icon)
                node['children'].append(key)
                self.__data_component[key] = {'kind': kind, 'path': path, 'file': file, 'children': None}
          else:  # here we search for new files and new folders that were never expanded
            node['children'] = []
            for item in files:
              kind, path, file = item.is_dir(), str(item.parent), item.name
              if kind == FILE and str(item) in displayed_files:
                continue
              key = new_key()
              displayed_files.add(str(item))
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

    capture.disable_proxy()
    self.__window.close()
