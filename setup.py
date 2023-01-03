#!/usr/bin/env python

from distutils.core import setup

try:
  from flat import __version__
except:
  print('Could not import package version')
  __version__ = "0.0.0"

setup(
  name='flat',
  version=__version__,
  description='PyMol Script Library/Collection/Repository',
  author='',
  author_email='',
  url='',
  packages=['flat'],
  package_dir='',
  package_data='',
  requires=['pymol'],
)