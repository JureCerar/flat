# Copyright (C) 2023-2024 Jure Cerar
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""
:mod:`flat.plotting`
====================
Module for plotting data within PyMOL.

.. note::

    Module has optional dependency to `mplcursors` python package.
    Install package with `pip install mplcursors`.
"""

from pymol import cmd, CmdException
import numpy as np 
import warnings
from pathlib import Path


def _showfigure(fig, filename, quiet):
    """Helper function for plot commands"""
    if not quiet:
        fig.show()
    if filename:
        fig.savefig(filename)
        if not quiet:
            print(f"Plot written to: {filename!r}")


def _get_model_color(selection, *, _self=cmd):
    """Get model color as RGB hex string to be used with matplotlib."""
    colors = []
    for guide in ("guide", "elem C", "all"):
        _self.iterate(
            f"({selection}) & {guide}",
            "colors.append(color)",
            space=locals()
        )
        if colors:
            break
    else:
        colors = ["gray",]

    # Find most frequent color
    color = max((colors.count(color), color) for color in colors)[1]

    if color >= 0x40000000:
        color = "0x%06x" % (color & 0xFFFFFF)
    
    # Return RGB hex
    color = _self.get_color_tuple(color)
    return "#%02x%02x%02x" % tuple(int(0xFF * v) for v in color)


@cmd.extend
def plot(expression="b", selection="all", fmt=None, byres=True,
         filename=None, *, quiet=1, _self=cmd, **kwargs):
    """
    DESCRIPTION
        Plot selected property or expression with matplotlib
    USAGE
        plot [ expression [, selection [, fmt, [, byres, [, filename ]]]]]
    ARGUMENTS
        expression : str, default = 'b'
            Atom property or expression to plot.
        selection : str, default = 'all'
            Atom selection.
        fmt : str, default = None
            Plotting format (see matplotlib documentation). By default plot as bars. 
        byres : bool, default = True
            Plot properties per-residue instead of per-atom.
        filename : str, optional
            Save figure to file.
    RETURNS
        : matplotlib.figure.Figure
            Figure object.
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    handles = list()
    labels = list()

    xlim = (np.inf, -np.inf)
    ylim = (np.inf, -np.inf)

    def _index_to_hex(index):
        """From PyMOL color index to hex"""
        val = _self.get_color_tuple(index)
        return "#%02x%02x%02x" % tuple(int(0xFF * v) for v in val)

    for model in _self.get_object_list(selection):
        for chain in _self.get_chains(model):
            colors = list()
            x_list = list()
            y_list = list()
            label_list = list()

            if int(byres):
                _self.iterate(
                    f"guide & ({selection}) & /{model}//{chain}",
                    f"x_list.append(int(resi));" +
                    f"y_list.append({expression});" +
                    "colors.append(color);" +
                    "label_list.append(f'/{model}/{segi}/{chain}/{resn}`{resi}');",
                    space=locals(),
                )
            else:
                _self.iterate(
                    f"{selection} & /{model}//{chain}",
                    f"x_list.append(index);" +
                    f"y_list.append({expression});" +
                    "colors.append(color);" +
                    "label_list.append(f'/{model}/{segi}/{chain}/{resn}`{resi}/{name}');",
                    space=locals(),
                )

            colors = [_index_to_hex(c) for c in colors]

            if not fmt:
                handle = ax.bar(x_list, y_list, width=1, color=colors,
                                label=f"{model}/{chain}", **kwargs)
            else:
                # Find most prominent color
                color = max((colors.count(c), c) for c in colors)[1]
                handle, = ax.plot(x_list, y_list, fmt, color=color,
                                  label=f"{model}/{chain}", **kwargs)
                
            xlim = min(xlim[0], min(x_list)), max(xlim[1], max(x_list))
            ylim = min(ylim[0], min(y_list)), max(ylim[1], max(y_list))

            handles.append(handle)
            labels.append(label_list)

    labels_dict = dict(zip(handles, labels))

    try:
        import mplcursors

        cursor = mplcursors.cursor(ax, hover=False)
        @cursor.connect("add")
        def on_click(sel):
            # BUG: Pymol return sel.index as numpy.float64?
            sele = labels_dict[sel.artist][int(sel.index)]
            sel.annotation.set_text(sele)
            sel.annotation.get_bbox_patch().set(fc="white", alpha=0.8)
            _self.select(sele)
            # _self.zoom(sele)

    except ImportError as error:
        warnings.warn(error)

    # Set custom limits for x- and y-axis
    yrange = ylim[1] - ylim[0]

    ax.set_xlim(xlim[0] - 1, xlim[1] + 1)
    ax.set_ylim(ylim[0] - 0.05 * yrange, ylim[1] + 0.05 * yrange)

    # Plot grid and title
    # ax.set_xlabel("Index")
    # ax.set_ylabel(expression)
    ax.grid()
    fig.legend()
    _showfigure(fig, filename, quiet)

    return fig


@cmd.extend
def plot_contacts(selection="guide", metric="euclidean", *,
                  state=-1, filename=None, quiet=1, _self=cmd):
    """
    DESCRIPTION
        Plot a contact map.
    USAGE
        plot_contacts [ selection [, metric [, state [, filename ]]]]
    ARGUMENTS
        selection : str, default = 'guide'
            Atom selection.
        metric : str, default = 'euclidean'
            Metric for distance matrix.
        state : int, default = -1
            State index or all states if state=0 {default: -1}
        filename : str, optional
            Save figure to file.
    RETURNS
        : matplotlib.figure.Figure
            Figure object.
    SOURCE
        From PSICO (c) 2011-2012 Thomas Holder, MPI for Developmental Biology
    """
    import matplotlib.pyplot as plt
    from .numeric import pdist_squareform

    X = _self.get_coords(selection, int(state))
    if X is None:
        raise CmdException("No coordinates in selection")

    dist_mat = pdist_squareform(X, metric)

    fig, ax = plt.subplots()
    mappable = ax.pcolormesh(dist_mat, cmap="viridis_r")
    fig.colorbar(mappable, ax=ax)

    _showfigure(fig, filename, quiet)

    return fig


@cmd.extend
def plot_ramachandran(selection="guide", fmt=".", state=-1, ref=1, 
                      filename=None, *, quiet=1, _self=cmd, **kwargs):
    """
    DESCRIPTION
        Plot a Ramachandra dihedral map.
    USAGE
        plot_ramachandran [ selection [, fmt [, state [, ref [, filename ]]]]]
    ARGUMENTS
        selection : str, default = 'guide'
            Atom selection.
        fmt : str, default = '.'
            Plotting format (see matplotlib documentation).
        state : int, default = -1
            State index or all states if state=0 {default: -1}
        ref : bool, default = True
            Plot reference Ramachandra dihedrals.
        filename : str, optional
            Save figure to file.
    RETURNS
        : matplotlib.figure.Figure
            Figure object.
    SOURCE
        Reference Ramachandra plot was taken from MDAnalysis dihedral package:
        Michaud-Agrawal, J. Comput. Chem., 2011, doi:10.1002/jcc.21787
    SEE ALSO
        :func:`phi_psi`
    """
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()
    handles = list()
    labels = list()

    for model in _self.get_object_list(selection):
        for chain in _self.get_chains(model):
            color = _get_model_color(f"({selection}) & /{model}//{chain}",
                                     _self=_self)
            r = _self.get_phipsi(f"({selection}) & /{model}//{chain}",
                                 state, _self=_self)
            
            if not r or not isinstance(r, dict):
                warnings.warn(f"No dihderal angles found for: '{model}/{chain}'")
                continue

            residues = []
            for model, index in sorted(r):
                _self.iterate(
                    f"{model}`{index}",
                    "residues.append(f'/{model}/{segi}/{chain}/{resn}`{resi}')",
                    space=locals(),
                )

            xy = np.array(list(r.values())).T
            handle, = ax.plot(xy[0], xy[1], fmt, color=color,
                            label=f"{model}/{chain}", zorder=2, **kwargs)
            
            handles.append(handle)
            labels.append(residues)
            
    labels_dict = dict(zip(handles, labels))

    try:
        # Add interactive cursors
        import mplcursors
        
        cursor = mplcursors.cursor(ax, hover=True)
        @cursor.connect("add")
        def on_hover(sel):
            # BUG: Pymol return sel.index as numpy.float64?
            sele = labels_dict[sel.artist][int(sel.index)]
            sel.annotation.set_text(sele)
            sel.annotation.get_bbox_patch().set(fc="white", alpha=0.8)

    except ImportError as error:
        warnings.warn(error)

    if int(ref):
        # Load reference Ramachandra plot
        path = Path(__file__).parent.resolve()
        rama_ref = Path.joinpath(path, "data", "rama_ref.npy")
        levels = [1, 17, 15000]
        colors = ["#eeeeee", "#dfdfdf"]
        X, Y = np.meshgrid(np.arange(-180, 180, 4), np.arange(-180, 180, 4))
        Z = np.load(rama_ref)
        ax.contourf(X, Y, Z, levels=levels, colors=colors, zorder=0)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-180, 181, 60))
    ax.set_xlabel(r"$\phi$ [deg]")
    ax.set_ylabel(r"$\psi$ [deg]")

    ax.grid(zorder=1)
    fig.legend()
    _showfigure(fig, filename, quiet)

    return fig


# Autocompletion
cmd.auto_arg[0].update({
    "plot": cmd.auto_arg[0]["spectrum"],
    "plot_contacts": cmd.auto_arg[0]["zoom"],
    "plot_ramachandran": cmd.auto_arg[0]["zoom"],
})
cmd.auto_arg[1].update({
    "plot": cmd.auto_arg[0]["zoom"],
})


