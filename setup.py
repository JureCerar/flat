from setuptools import find_packages, setup

try:
  from flat import __version__
except:
  print('Could not import package version')
  __version__ = "0.0.0"

setup(
  name='flat',
  packages=find_packages(),
  version=__version__,
  license='GNU GPL 3.0',
  description='PyMol Script Library/Collection/Repository',
  author='Jure Cerar',
  url='',
  requires=['pymol'],
)
