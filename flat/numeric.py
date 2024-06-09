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

def pdist_squareform(X, metric="euclidean"):
    """Distance matrix. Uses scipy if available (faster) or numpy."""
    try:
        from scipy.spatial.distance import pdist, squareform
        return squareform(pdist(X, metric))
    
    except ModuleNotFoundError:
        import numpy

        X = numpy.asarray(X)
        Q = (X**2).sum(1)
        S = numpy.add.outer(Q, Q)
        D = (S - 2 * (X @ X.T)).clip(0)  # sqeuclidean

        if metric == "sqeuclidean":
            return D
        elif metric == "euclidean":
            return D**0.5

        raise ValueError(f"Unknown Distance Metric: {metric}")
