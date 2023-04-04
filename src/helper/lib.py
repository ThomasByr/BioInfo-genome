import sys
import os
import inspect
import threading
from typing import Any, NoReturn

from dotenv import load_dotenv
from termcolor import colored

__all__ = ['debug', 'info', 'error', 'panic']

load_dotenv()

__debug = os.getenv('DEBUG')
__is_debug: bool = False
if __debug is not None:
  if isinstance(__debug, str):
    if __debug.lower() in {'true', '1'}:
      __is_debug = True
  elif isinstance(__debug, int):
    if __debug == 1:
      __is_debug = True
  elif isinstance(__debug, bool):
    __is_debug = __debug

__print_msg_lock = threading.Lock()


def __print_msg(msg: str | Any) -> None:
  __print_msg_lock.acquire()
  print(msg if msg is not None else '', file=sys.stderr, flush=True)
  __print_msg_lock.release()


def debug(msg: str | Any = None) -> None:
  """
  print debug message to stderr
  
  ## Parameters
  ```py
  >>> msg : str
  ```
  string to print
  """
  if not __is_debug:
    return
  __print_msg(colored('  [debug] ', 'green') + msg)


def info(msg: str | Any = None) -> None:
  """
  print info message to stderr
  
  ## Parameters
  ```py
  >>> msg : str
  ```
  string to print
  """
  __print_msg(colored('   [info] ', 'blue') + msg)


def error(msg: str | Any = None) -> None:
  """
  print error message to stderr
  
  ## Parameters
  ```py
  >>> msg : str
  ```
  string to print
  """
  __print_msg(colored('  [error] ', 'yellow') + msg)


def panic(msg: str | Any = None) -> NoReturn:
  """
  print error message to stderr and exit\\
  this function does not return

  ## Parameters
  ```py
  >>> msg : str
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
