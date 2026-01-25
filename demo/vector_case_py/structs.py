# -*- coding: utf-8 -*-
"""
Data structures for VECTOR-CASE calculations.

This module defines the core data structures used for mutual inductance
calculations, including coil geometry and precision parameters.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Coil:
    """
    Representation of a coaxial coil for inductance calculations.

    A coil is defined by its geometry (inner radius, outer radius, length)
    and the number of wire turns.

    Parameters
    ----------
    inner_radius : float
        Inner radius of the coil in meters. Must be positive.
    outer_radius : float
        Outer radius of the coil in meters. Must be greater than ``inner_radius``.
    length : float
        Axial length of the coil in meters. Must be positive.
    num_turns : float
        Number of wire turns. Must be positive.

    Raises
    ------
    ValueError
        If any geometric constraint is violated.

    Examples
    --------
    >>> coil = Coil(inner_radius=0.1, outer_radius=0.2, length=0.1, num_turns=100)
    >>> coil.inner_radius
    0.1
    """

    inner_radius: float
    outer_radius: float
    length: float
    num_turns: float

    def __post_init__(self) -> None:
        """Validate coil parameters after initialization."""
        if self.inner_radius <= 0:
            raise ValueError(f"inner_radius must be positive, got {self.inner_radius}")
        if self.outer_radius <= 0:
            raise ValueError(f"outer_radius must be positive, got {self.outer_radius}")
        if self.inner_radius >= self.outer_radius:
            raise ValueError(
                f"inner_radius ({self.inner_radius}) must be less than "
                f"outer_radius ({self.outer_radius})"
            )
        if self.length <= 0:
            raise ValueError(f"length must be positive, got {self.length}")
        if self.num_turns <= 0:
            raise ValueError(f"num_turns must be positive, got {self.num_turns}")


@dataclass
class SumPrecision:
    """
    Precision parameters for the series expansion computation.

    Controls the number of terms used in the series expansion for mutual
    inductance calculation. Higher values yield more accurate results at
    the cost of increased computation time.

    Parameters
    ----------
    k_terms : int
        Number of terms for the k-index series. Default is 32.
    l_terms : int
        Number of terms for the l-index series. Default is 32.
    n_terms : int
        Number of terms for the n-index series. Default is 64.

    Notes
    -----
    For optimal accuracy:

    - ``k_terms == l_terms`` is recommended in general.
    - ``n_terms == 2 * k_terms`` is recommended for near-field calculations.
    - ``n_terms == k_terms`` is recommended for far-field calculations.

    Examples
    --------
    >>> precision = SumPrecision(k_terms=32, l_terms=32, n_terms=64)
    >>> precision.k_terms
    32
    """

    k_terms: int = 32
    l_terms: int = 32
    n_terms: int = 64

    def __post_init__(self) -> None:
        """Validate precision parameters after initialization."""
        if self.k_terms <= 0:
            raise ValueError(f"k_terms must be positive, got {self.k_terms}")
        if self.l_terms <= 0:
            raise ValueError(f"l_terms must be positive, got {self.l_terms}")
        if self.n_terms <= 0:
            raise ValueError(f"n_terms must be positive, got {self.n_terms}")
