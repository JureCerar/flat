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
