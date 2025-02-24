"""
ls --- :mod:`lipid_order.analysis.ls`
===========================================================

This module contains the :class:`ls` class.

"""
from typing import Union, TYPE_CHECKING

from MDAnalysis.analysis.base import AnalysisBase
import numpy as np

if TYPE_CHECKING:
    from MDAnalysis.core.universe import Universe, AtomGroup


class ls(AnalysisBase):
    """ls class.

    This class is used to perform analysis on a trajectory.

    Parameters
    ----------
    universe_or_atomgroup: :class:`~MDAnalysis.core.universe.Universe` or :class:`~MDAnalysis.core.groups.AtomGroup`
        Universe or group of atoms to apply this analysis to.
        If a trajectory is associated with the atoms,
        then the computation iterates over the trajectory.
    select: str
        Selection string for atoms to extract from the input Universe or
        AtomGroup

    Attributes
    ----------
    universe: :class:`~MDAnalysis.core.universe.Universe`
        The universe to which this analysis is applied
    atomgroup: :class:`~MDAnalysis.core.groups.AtomGroup`
        The atoms to which this analysis is applied
    results: :class:`~MDAnalysis.analysis.base.Results`
        results of calculation are stored here, after calling
        :meth:`ls.run`
    start: Optional[int]
        The first frame of the trajectory used to compute the analysis
    stop: Optional[int]
        The frame to stop at for the analysis
    step: Optional[int]
        Number of frames to skip between each analyzed frame
    n_frames: int
        Number of frames analysed in the trajectory
    times: numpy.ndarray
        array of Timestep times. Only exists after calling
        :meth:`ls.run`
    frames: numpy.ndarray
        array of Timestep frame indices. Only exists after calling
        :meth:`ls.run`
    """

    def __init__(
        self,
        universe_or_atomgroup: Union["Universe", "AtomGroup"],
        select: str = "all",
        # TODO: add your own parameters here
        **kwargs
    ):
        # the below line must be kept to initialize the AnalysisBase class!
        super().__init__(universe_or_atomgroup.trajectory, **kwargs)
        # after this you will be able to access `self.results`
        # `self.results` is a dictionary-like object
        # that can should used to store and retrieve results
        # See more at the MDAnalysis documentation:
        # https://docs.mdanalysis.org/stable/documentation_pages/analysis/base.html?highlight=results#MDAnalysis.analysis.base.Results

        self.universe = universe_or_atomgroup.universe
        self.atomgroup = universe_or_atomgroup.select_atoms(select)

    def _prepare(self):
        """Set things up before the analysis loop begins"""
        # This is an optional method that runs before
        # _single_frame loops over the trajectory.
        # It is useful for setting up results arrays
        # For example, below we create an array to store
        # the number of atoms with negative coordinates
        # in each frame.
        self.results.is_negative = np.zeros(
            (self.n_frames, self.atomgroup.n_atoms),
            dtype=bool,
        )

    def _single_frame(self):
        """Calculate data from a single frame of trajectory"""
        # This runs once for each frame of the trajectory
        # It can contain the main analysis method, or just collect data
        # so that analysis can be done over the aggregate data
        # in _conclude.

        # The trajectory positions update automatically
        negative = self.atomgroup.positions < 0
        # You can access the frame number using self._frame_index
        self.results.is_negative[self._frame_index] = negative.any(axis=1)

    def _conclude(self):
        """Calculate the final results of the analysis"""
        # This is an optional method that runs after
        # _single_frame loops over the trajectory.
        # It is useful for calculating the final results
        # of the analysis.
        # For example, below we determine the
        # which atoms always have negative coordinates.
        self.results.always_negative = self.results.is_negative.all(axis=0)
        always_negative_atoms = self.atomgroup[self.results.always_negative]
        self.results.always_negative_atoms = always_negative_atoms
        self.results.always_negative_atom_names = always_negative_atoms.names

        # results don't have to be arrays -- they can be any value, e.g. floats
        self.results.n_negative_atoms = self.results.is_negative.sum(axis=1)
        self.results.mean_negative_atoms = self.results.n_negative_atoms.mean()
