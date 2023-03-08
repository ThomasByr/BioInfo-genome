from src import Tree

from src.helper import info, debug

if __name__ == '__main__':
  overview = Tree()
  eukaryotes = Tree('eukaryotes')
  plasmids = Tree('plasmids')
  prokaryotes = Tree('prokaryotes')
  viruses = Tree('viruses')

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

  main_tree.save('main_tree')
  debug(f'wrote to file: {main_tree.filepath}')
