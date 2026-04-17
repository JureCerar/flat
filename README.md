<p align="center">
    <img src="docs/static/flatcat.png" alt="adjust" width="300">
</p> 

# FlatCAT - PyMOL Script Collection

**FlatCAT** is personal collection of scripts and utilities that extend functionality of PyMOL molecular visualization tool. The inspiration behind the library is to learn how to work with PyMOL and create a collection of useful functions for easy and convenient use. It is heavily inspired and borrows (i.e. *shamelessly steals*) ideas from [Pymol ScrIpt COllection (PSICO)](https://github.com/speleo3/pymol-psico).

## Installation

The package can be installed by cloning the repository and running:

```bash
pip install .
```

## Usage

Once installed, you can import and initialize package in PyMOL by typing:

```python
import flat.fullinit
```

> [!NOTE]
> To import **FlatCAT** package in PyMOL at startup, add this to your `~/.pymolrc.py` file or add it in PyMOL under `File` → `Edit pymolrc`.

You can then use built-in help in PyMOL to get command documentation:

```text
PyMOL> help flat
```
## Documentation

Documentation is available online at [Read the Docs](https://flatcat.readthedocs.io/en/latest/).

## Getting help

Please report bugs and feature requests for **FlatCAT** through the [Issue Tracker](https://github.com/JureCerar/flat/issues).

## Contributing

All contributions are welcome. To contribute code, submit a pull request against the master branch in the repository.

## License

This program is licensed under the GNU General Public License v3.0

Copyright (C) 2023-2026 [Jure Cerar](https://github.com/JureCerar)

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.