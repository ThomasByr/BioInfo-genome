import os

from src import Tree, GenomeGUI

from src.helper import info, debug

if __name__ == '__main__':
  overview = Tree()

  # Tree.update_ids()
  gui = GenomeGUI('src')
  gui.run()
  # Tree.clean_folders()
