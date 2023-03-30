import os

from pathlib import Path
import PySimpleGUI as sg

from ..helper import panic, error, info, debug

__all__ = ['GenomeGUI']

folder_icon = b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAACXBIWXMAAAsSAAALEgHS3X78AAABnUlEQVQ4y8WSv2rUQRSFv7vZgJFFsQg2EkWb4AvEJ8hqKVilSmFn3iNvIAp21oIW9haihBRKiqwElMVsIJjNrprsOr/5dyzml3UhEQIWHhjmcpn7zblw4B9lJ8Xag9mlmQb3AJzX3tOX8Tngzg349q7t5xcfzpKGhOFHnjx+9qLTzW8wsmFTL2Gzk7Y2O/k9kCbtwUZbV+Zvo8Md3PALrjoiqsKSR9ljpAJpwOsNtlfXfRvoNU8Arr/NsVo0ry5z4dZN5hoGqEzYDChBOoKwS/vSq0XW3y5NAI/uN1cvLqzQur4MCpBGEEd1PQDfQ74HYR+LfeQOAOYAmgAmbly+dgfid5CHPIKqC74L8RDyGPIYy7+QQjFWa7ICsQ8SpB/IfcJSDVMAJUwJkYDMNOEPIBxA/gnuMyYPijXAI3lMse7FGnIKsIuqrxgRSeXOoYZUCI8pIKW/OHA7kD2YYcpAKgM5ABXk4qSsdJaDOMCsgTIYAlL5TQFTyUIZDmev0N/bnwqnylEBQS45UKnHx/lUlFvA3fo+jwR8ALb47/oNma38cuqiJ9AAAAAASUVORK5CYII='
file_icon = b'iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAYAAAAf8/9hAAAACXBIWXMAAAsSAAALEgHS3X78AAABU0lEQVQ4y52TzStEURiHn/ecc6XG54JSdlMkNhYWsiILS0lsJaUsLW2Mv8CfIDtr2VtbY4GUEvmIZnKbZsY977Uwt2HcyW1+dTZvt6fn9557BGB+aaNQKBR2ifkbgWR+cX13ubO1svz++niVTA1ArDHDg91UahHFsMxbKWycYsjze4muTsP64vT43v7hSf/A0FgdjQPQWAmco68nB+T+SFSqNUQgcIbN1bn8Z3RwvL22MAvcu8TACFgrpMVZ4aUYcn77BMDkxGgemAGOHIBXxRjBWZMKoCPA2h6qEUSRR2MF6GxUUMUaIUgBCNTnAcm3H2G5YQfgvccYIXAtDH7FoKq/AaqKlbrBj2trFVXfBPAea4SOIIsBeN9kkCwxsNkAqRWy7+B7Z00G3xVc2wZeMSI4S7sVYkSk5Z/4PyBWROqvox3A28PN2cjUwinQC9QyckKALxj4kv2auK0xAAAAAElFTkSuQmCC'


class GenomeGUI:

  def __init__(self, results_folder: str = 'Results'):
    self.__window: sg.Window = None
    self.__results_folder = results_folder
    self.__font = ('Helvetica', 11)
    self.__monospace_font = ('Courier', 7)
    self.__file_tree: list[list[sg.Tree]] = None
    self.__tree_data: sg.TreeData = None
    self.__selected_organism: str = None
    self.__info = [[sg.Text('Information', font=self.__font)], [sg.Multiline(
        size=(50, 10),
        key='-INFO-',
    )]]
    self.__log = [[sg.Text('Logs', font=self.__font)],
                  [
                      sg.Multiline(
                          key='-LOG-',
                          size=(60, 20),
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
                  sg.Button('Cancel', font=self.__font, size=(10, 2)),
                  sg.Button('Confirm Selection', font=self.__font, size=(10, 2)),
              ]]
    return layout

  def __build_file_tree(self):
    self.__file_tree = [[
        sg.Tree(
            data=self.__build_tree(),
            headings=[],
            auto_size_columns=False,
            num_rows=20,
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
    sg.theme('DarkTeal9')
    self.__build_tree()
    self.__build_file_tree()
    self.__build_window()

    def new_key() -> int:
      key = 1
      while key in data:
        key += 1
      return key

    tree = self.__window['-TREE-']
    if isinstance(tree, sg.ErrorElement):
      panic('Tree not found')

    tree.Widget.configure(show='tree')
    tree.bind('<Double-1>', "DOUBLE-CLICK-")
    tree.bind('<Button-3>', "RIGHT-CLICK-")

    DIR, FILE = True, False  #pylint: disable=unused-variable
    data = {
        0: {
            'kind': DIR,
            'path': '',
            'file': os.path.basename(self.__results_folder),
            'children': None
        },
    }

    while True:
      event, values = self.__window.read()
      if event in (sg.WIN_CLOSED, 'Exit'):
        break

      if event == 'Confirm Selection':
        try:
          selected = values['-TREE-'][0]
          self.__selected_organism = data[selected]['file']
          info(f'Selected organism: {self.__selected_organism}')
        except IndexError:
          error('Nothing selected in the tree')

      if event == '-TREE-DOUBLE-CLICK-':
        parent_key = values['-TREE-'][0]
        node = data[parent_key]

        if node['kind'] == DIR and node['children'] is None:
          parent_path = Path(node['path']).joinpath(node['file'])
          files: list
          try:
            files = sorted(list(parent_path.iterdir()), key=lambda file: file.is_file())
          except PermissionError:
            error(f'Permission denied: {parent_path}')
            continue

          node['children'] = []
          for item in files:
            key = new_key()
            kind, path, file = item.is_dir(), str(item.parent), item.name
            self.__tree_data.insert(parent_key,
                                    key,
                                    str(file), [],
                                    icon=folder_icon if kind == DIR else file_icon)
            node['children'].append(key)
            data[key] = {'kind': kind, 'path': path, 'file': file, 'children': None}

          debug(f'Loaded {len(files)} items from {parent_path}')
          tree.update(values=self.__tree_data)
          iid = tree.KeyToID[parent_key]
          tree.Widget.see(iid)

    self.__window.close()
