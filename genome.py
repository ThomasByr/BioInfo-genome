
from src import Tree, GenomeGUI

if __name__ == '__main__':
  overview = Tree()
  overview.build()

  # Tree.update_ids()
  gui = GenomeGUI()
  gui.run()
  # Tree.clean_folders()
