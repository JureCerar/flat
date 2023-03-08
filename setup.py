<<<<<<< HEAD
from setuptools import find_packages, setup

setup(
  name='flat',
  packages=find_packages(),
  version='0.0.0',
  description='',
  author='Jure Cerar',
  license='GNU GPL 3.0',
  install_requires=[],
  setup_requires=['pytest-runner'],
  tests_require=['pytest'],
  test_suite='tests',
)
=======
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
>>>>>>> 88f069c7cace1fd1d18c6d01a83bba514d9292da
