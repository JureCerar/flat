from pymol import cmd, menu

# @cmd.extend
# def cbc(color="green"):
#     """ Color by chain """
#     pass

# @cmd.extend
# def cbe(color="green"):
#     """ Color by element """
#     pass

# @cmd.extend
# def cbattr(color="green"):
#     """ Color by residue attribute """
#     pass


class _Palette:
    """ Dummy class """
    pass


class _Colors():
    """ Class for storing color palettes. PyMOL API only. """

    def __init__(self, *args):
        self.colors = list(args)
        self.index = 0
        return

    def __repr__(self):
        return " ".join(self.colors)

    def replace(self, *args):
        return " ".join(self.colors)

    def __iter__(self):
        return self

    def __next__(self):
        try:
            result = self.colors[self.index]
        except IndexError:
            raise StopIteration
        self.index += 1
        return result

    def list(self):
        return self.colors


# Color palettes. PyMOL API only
palette = _Palette()

# Novartis colors
palette.novartis = _Colors("novartis", "carmine", "sienna", "apricot")
palette.novartis_rev = _Colors("apricot", "sienna", "carmine", "novartis")
cmd.set_color("novartis", [4, 96, 169])
cmd.set_color("carmine", [141, 21, 27])
cmd.set_color("sienna", [231, 74, 33])
cmd.set_color("apricot", [236, 154, 30])
menu.all_colors_list.append(
    ('novartis', [
        ('036', 'novartis'),
        ('501', 'carmine'),
        ('921', 'sienna'),
        ('961', 'apricot'),
    ]),
)

# AlphaFold colors
palette.alphafold = _Colors("af_high", "af_medium", "af_low", "af_crit")
palette.alphafold_rev = _Colors("af_crit", "af_low", "af_medium", "af_high")
cmd.set_color("af_high", [0, 83, 214])
cmd.set_color("af_medium", [101, 203, 243])
cmd.set_color("af_low", [255, 219, 19])
cmd.set_color("af_crit", [255, 125, 69])
menu.all_colors_list.append(
    ('alphafold', [
        ('038', 'af_high'),
        ('379', 'af_medium'),
        ('980', 'af_low'),
        ('942', 'af_crit'),
    ]),
)
