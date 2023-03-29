import os

from src import Tree, GenomeGUI

from src.helper import info, debug

if __name__ == '__main__':
  overview = Tree()
  eukaryotes = Tree('eukaryotes')
  plasmids = Tree('plasmids')
  prokaryotes = Tree('prokaryotes')
  viruses = Tree('viruses')

  # if "main_tree.json" exists in data directory, load it
  # else, build the tree from scratch
  should_build = False
  main_tree: Tree = None
  try:
    os.open('data/main_tree.json', os.O_RDONLY)
    os.open('data/overview.json', os.O_RDONLY)
    main_tree = Tree.load('main_tree')
    overview = Tree.load('overview')
    info(f'overview loaded from file: {overview.filepath}')
    info(f'main tree loaded from file: {main_tree.filepath}')
  except FileNotFoundError:
    info('main tree not found, building from scratch')
    should_build = True

  if should_build:
    overview.build()
    info(f'overview built with {(n:=len(overview))} entr{"y" if n <= 1 else "ies"}')
    eukaryotes.build()
    info(f'eukaryotes built with {(n:=len(eukaryotes))} entr{"y" if n <= 1 else "ies"}')
    plasmids.build()
    info(f'plasmids built with {(n:=len(plasmids))} entr{"y" if n <= 1 else "ies"}')
    prokaryotes.build()
    info(f'prokaryotes built with {(n:=len(prokaryotes))} entr{"y" if n <= 1 else "ies"}')
    viruses.build()
    info(f'viruses built with {(n:=len(viruses))} entr{"y" if n <= 1 else "ies"}')
    info()

    main_tree = eukaryotes.union(plasmids).union(prokaryotes).union(viruses)
    info(f'main tree built with {(n := len(main_tree))} kingdom{"" if n <= 1 else "s"}')
    info(f'kingdoms: {", ".join(main_tree.tree.keys())}')

    overview.save('overview')
    debug(f'wrote to file: {overview.filepath}')
    main_tree.save('main_tree')
    debug(f'wrote to file: {main_tree.filepath}')

  overview.to_folders()
  debug(f'overview system folders created')
  # Tree.update_ids()
  gui = GenomeGUI()
  gui.run()
  # Tree.clean_folders()
