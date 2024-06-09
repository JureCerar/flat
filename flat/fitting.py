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

from pymol import cmd, CmdException

class MatchMaker(object):
    """
    DESCRIPTION
        API only. Matches two atom selections and provides two matched
        subselections with equal atom count. May involve temporary objects
        or named selections which will be automatically deleted.
    ARGUMENTS
        mobile = string: first atom selection
        target = string: second atom selection
        match = string: method how to match atoms
            * none: (dummy)
            * in: match atoms by "in" operator
            * like: match atoms by "like" operator
            * align: match atoms by cmd.align (without refinement)
            * super: match atoms by cmd.super (without refinement)
            * <name of alignment object>: use given alignment
    RESULT
        Properties "mobile" and "target" hold the matched subselections as
        selection strings.
    SOURCE
        From PSICO (c) 2011 Thomas Holder, MPI for Developmental Biology
    """

    def __init__(self, mobile, target, match, *, autodelete=True, _self=cmd):
        self._self = _self
        self.autodelete = autodelete
        self.temporary = []

        if match == "none":
            self.mobile = mobile
            self.target = target
        elif match in ["in", "like"]:
            self.mobile = "(%s) %s (%s)" % (mobile, match, target)
            self.target = "(%s) %s (%s)" % (target, match, mobile)
        elif match in ["align", "super"]:
            self.align(mobile, target, match)
        elif match in _self.get_names("all") and _self.get_type(match) in ("object:", "object:alignment"):
            self.from_alignment(mobile, target, match)
        else:
            raise CmdException("unkown match method", match)

    def check(self):
        return self._self.count_atoms(self.mobile) == self._self.count_atoms(self.target)

    def align(self, mobile, target, match):
        """
        Align mobile to target using the alignment method given by "match"
        """
        aln_obj = self._self.get_unused_name("_")
        self.temporary.append(aln_obj)

        align = cmd.keyword[match][0]
        align(mobile, target, cycles=0, transform=0,
              object=aln_obj, _self=self._self)
        self._self.disable(aln_obj)

        self.from_alignment(mobile, target, aln_obj)

    def from_alignment(self, mobile, target, aln_obj):
        """
        Use alignment given by "aln_obj" (name of alignment object)
        """
        from .selecting import wait_for
        wait_for(aln_obj, _self=self._self)

        self.mobile = "(%s) and %s" % (mobile, aln_obj)
        self.target = "(%s) and %s" % (target, aln_obj)
        if self.check():
            return

        # difficult: if selections spans only part of the alignment or
        # if alignment object covers more than the two objects, then we
        # need to pick those columns that have no gap in any of the two
        # given selections

        mobileidx = set(self._self.index(mobile))
        targetidx = set(self._self.index(target))
        mobileidxsel = []
        targetidxsel = []

        for column in self._self.get_raw_alignment(aln_obj):
            mobiles = mobileidx.intersection(column)
            if len(mobiles) == 1:
                targets = targetidx.intersection(column)
                if len(targets) == 1:
                    mobileidxsel.extend(mobiles)
                    targetidxsel.extend(targets)

        self.mobile = self._self.get_unused_name("_mobile")
        self.target = self._self.get_unused_name("_target")
        self.temporary.append(self.mobile)
        self.temporary.append(self.target)

        mobile_objects = set(idx[0] for idx in mobileidxsel)
        target_objects = set(idx[0] for idx in targetidxsel)

        if len(mobile_objects) == len(target_objects) == 1:
            mobile_index_list = [idx[1] for idx in mobileidxsel]
            target_index_list = [idx[1] for idx in targetidxsel]
            self._self.select_list(
                self.mobile, mobile_objects.pop(), mobile_index_list, mode="index")
            self._self.select_list(
                self.target, target_objects.pop(), target_index_list, mode="index")
        else:
            self._self.select(self.mobile, " ".join("%s`%d" %
                              idx for idx in mobileidxsel))
            self._self.select(self.target, " ".join("%s`%d" %
                              idx for idx in targetidxsel))

    def __enter__(self):
        return self

    def __exit__(self, type_, value, traceback):
        self._cleanup()
        self.autodelete = False

    def __del__(self):
        self._cleanup()

    def _cleanup(self):
        if not self.autodelete:
            return
        for name in self.temporary:
            self._self.delete(name)
