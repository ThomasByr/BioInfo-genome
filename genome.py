
from src import Tree, GenomeGUI

if __name__ == '__main__':
  overview = Tree()
  overview.build()

  gui = GenomeGUI(overview)
  gui.run()
