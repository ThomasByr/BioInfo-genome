import sys
import os
import inspect
from typing import Any, NoReturn

from termcolor import colored

__all__ = ['debug', 'info', 'error', 'panic']


def _print_msg(msg: str | Any) -> None:
  print(msg if msg is not None else '', file=sys.stderr)
  sys.stderr.flush()


def debug(msg: str | Any = None) -> None:
  """
  print debug message to stderr
  
  ## Parameters
  ```py
    msg : str
  ```
    string to print
  """
  print(colored('  [debug]', 'green'), file=sys.stderr, end=' ')
  _print_msg(msg)


def info(msg: str | Any = None) -> None:
  """
  print info message to stderr
  
  ## Parameters
  ```py
    msg : str
  ```
    string to print
  """
  print(colored('   [info]', 'blue'), file=sys.stderr, end=' ')
  _print_msg(msg)


def error(msg: str | Any = None) -> None:
  """
  print error message to stderr
  
  ## Parameters
  ```py
    msg : str
  ```
    string to print
  """
  print(colored('  [error]', 'yellow'), file=sys.stderr, end=' ')
  _print_msg(msg)


def panic(msg: str | Any = None) -> NoReturn:
  """
  print error message to stderr and exit\\
  this function does not return

  ## Parameters
  ```py
    msg : str
  ```
    string to print
  ## Raises
  ```py
    RuntimeError
  ```
  """
  # yapf: disable
  output = [
    colored('  [panic] ', 'red'), msg, '\n',
    colored('       -> ', 'red'),
    f'stack trace: {os.path.basename(inspect.stack()[1].filename)}:'
    f'{inspect.stack()[1].function}:'f'{inspect.stack()[1].lineno}',
    '\n', colored('       -> ', 'red'),
    f'exception: {sys.exc_info()[1]}',
  ]
  # yapf: enable

  raise RuntimeError('\n' + (''.join(output)))
