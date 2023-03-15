# the definition of the GenomeGUI class

import PySimpleGUI as sg
from ..core.tree import *

__all__ = ['GenomeGUI']


class GenomeGUI:

  def __init__(self, tree: Tree):
    self.__window = None
    self.__tree = tree
    self.__data = sg.TreeData()
    self.__font = ('Helvetica', 11)
    self.__file_tree = None
    self.__info = [[sg.Text('Information', font=self.__font)],
                   [sg.Multiline(size=(50, 10), key='-INFO-', disabled=True)]]
    self.__log = [[sg.Text('Logs', font=self.__font)],
                  [sg.Multiline(size=(50, 10), key='-LOG-', disabled=True)]]
    self.__right_side = [[sg.Column(self.__info, vertical_alignment='top')], [sg.HSeparator()],
                         [sg.Column(self.__log, vertical_alignment='top')]]

  def __build_tree(self):
    if self.__tree is None:
      return
    for kingdom in self.__tree.tree.items():
      kingdom = kingdom[1]
      self.__data.Insert('', kingdom.name, kingdom.name, values=['', ''])
      for group in kingdom.groups.values():
        self.__data.Insert(kingdom.name, group.name, group.name, values=['', ''])
        for subgroup in group.subgroups.values():
          self.__data.Insert(group.name, subgroup.name, subgroup.name, values=['', ''])

  def __build_layout(self):
    layout = [[sg.Text('GENOME', font=self.__font)], [sg.HSeparator()],
              [
                  sg.Column(self.__file_tree, justification='center'),
                  sg.VSeparator(),
                  sg.Column(self.__right_side)
              ], [sg.HSeparator()],
              [sg.Button('Run', font=self.__font),
               sg.Button('Cancel', font=self.__font)]]
    return layout

  def __build_file_tree(self):
    self.__file_tree = [[
        sg.Tree(data=self.__data,
                headings=['Size', 'Date Modified'],
                auto_size_columns=True,
                num_rows=20,
                col0_width=20,
                key='-TREE-')
    ]]

  def __build_window(self):
    self.__window = sg.Window('Genome', self.__build_layout(), size=(800, 600))

  def run(self):
    # Define the theme
    sg.theme('DarkTeal9')
    self.__build_tree()
    self.__build_file_tree()
    self.__build_window()
    while True:
      event, values = self.__window.read()
      if event == sg.WIN_CLOSED:
        break
    self.__window.close()
