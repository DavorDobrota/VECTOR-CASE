# -*- coding: utf-8 -*-
"""
Mutual inductance calculator for VECTOR-CASE.

This module provides the main calculator class for computing mutual inductance
between coaxial coils using the VECTOR-CASE C library.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Any, Optional

from .compiler import get_library
from .structs import Coil, SumPrecision

if TYPE_CHECKING:
    from types import ModuleType


class MutualInductanceCalculator:
    """
    Calculator for mutual inductance between coaxial coils.

    This class provides a high-level interface to the VECTOR-CASE C library
    for computing mutual inductance between two coaxial coils at a given
    separation distance.

    Parameters
    ----------
    precision : SumPrecision, optional
        Precision parameters for the series expansion. If None, uses defaults.
    auto_compile : bool, optional
        If True, automatically compile the library if needed. Default is True.

    Attributes
    ----------
    precision : SumPrecision
        The precision parameters used for calculations.

    Examples
    --------
    Basic calculation::

        from vector_case_py import Coil, MutualInductanceCalculator

        calculator = MutualInductanceCalculator()
        coil1 = Coil(inner_radius=0.1, outer_radius=0.2, length=0.1, num_turns=100)
        coil2 = Coil(inner_radius=0.3, outer_radius=0.4, length=0.1, num_turns=100)
        inductance = calculator.calculate(coil1, coil2, distance=0.1)

    With custom precision::

        from vector_case_py import SumPrecision

        precision = SumPrecision(k_terms=64, l_terms=64, n_terms=128)
        calculator = MutualInductanceCalculator(precision=precision)
    """

    def __init__(
        self,
        precision: Optional[SumPrecision] = None,
        auto_compile: bool = True,
    ) -> None:
        """
        Initialize the calculator with the given precision.

        Parameters
        ----------
        precision : SumPrecision, optional
            Precision parameters. Defaults to SumPrecision().
        auto_compile : bool, optional
            Whether to compile library on init. Default is True.
        """
        self.precision = precision or SumPrecision()
        self._lib: Optional["ModuleType"] = None

        if auto_compile:
            self._ensure_library()

    def _ensure_library(self) -> None:
        """Ensure the C library is loaded."""
        if self._lib is None:
            self._lib = get_library()

    @property
    def lib(self) -> "ModuleType":
        """
        Access the underlying C library module.

        Returns
        -------
        ModuleType
            The loaded vector_case module with ``ffi`` and ``lib`` attributes.
        """
        self._ensure_library()
        assert self._lib is not None
        return self._lib

    def calculate(
        self,
        coil1: Coil,
        coil2: Coil,
        distance: float,
        verbose: bool = False,
    ) -> float:
        """
        Calculate the mutual inductance between two coils.

        This method automatically selects the appropriate calculation method
        (near-field or far-field) based on the coil geometry and distance.

        Parameters
        ----------
        coil1 : Coil
            The first coil.
        coil2 : Coil
            The second coil.
        distance : float
            Axial distance between the coils in meters. Must be positive.
        verbose : bool, optional
            If True, print optimization details. Default is False.

        Returns
        -------
        float
            The mutual inductance in Henries.

        Raises
        ------
        ValueError
            If the distance is not positive.
        RuntimeError
            If the calculation fails (returns -1.0 from C library).

        Examples
        --------
        >>> calculator = MutualInductanceCalculator()  # doctest: +SKIP
        >>> coil1 = Coil(0.1, 0.2, 0.1, 100)
        >>> coil2 = Coil(0.3, 0.4, 0.1, 100)
        >>> inductance = calculator.calculate(coil1, coil2, distance=0.1)  # doctest: +SKIP
        """
        if distance <= 0:
            raise ValueError(f"distance must be positive, got {distance}")

        result = self.lib.lib.calculate_mutual_inductance_raw(
            coil1.inner_radius,
            coil1.outer_radius,
            coil1.length,
            coil1.num_turns,
            coil2.inner_radius,
            coil2.outer_radius,
            coil2.length,
            coil2.num_turns,
            distance,
            self.precision.k_terms,
            self.precision.l_terms,
            self.precision.n_terms,
            verbose,
        )

        if result == -1.0:
            raise RuntimeError(
                "Mutual inductance calculation failed. Check input parameters."
            )

        return result

    def calculate_near(
        self,
        coil1: Coil,
        coil2: Coil,
        distance: float,
        z_parameter: float,
    ) -> float:
        """
        Calculate mutual inductance using the near-field method.

        This method is optimized for coils that are close together. It requires
        an explicit Z parameter which affects the convergence of the series.

        Parameters
        ----------
        coil1 : Coil
            The first coil.
        coil2 : Coil
            The second coil.
        distance : float
            Axial distance between the coils in meters.
        z_parameter : float
            The Z parameter for the near-field calculation.

        Returns
        -------
        float
            The mutual inductance in Henries.

        Notes
        -----
        This is a lower-level method. For most use cases, prefer
        :meth:`calculate` which automatically selects the appropriate
        method and optimizes the Z parameter.

        See Also
        --------
        calculate : Automatic method selection.
        calculate_far : Far-field calculation method.
        """
        coil_data = self._create_coil_data(coil1, coil2)
        precision_data = self._create_precision_data()

        return self.lib.lib.calculate_mutual_inductance_near(
            coil_data[0], precision_data[0], distance, z_parameter
        )

    def calculate_far(
        self,
        coil1: Coil,
        coil2: Coil,
        distance: float,
    ) -> float:
        """
        Calculate mutual inductance using the far-field method.

        This method is optimized for coils that are far apart, where the
        far-field approximation converges quickly.

        Parameters
        ----------
        coil1 : Coil
            The first coil.
        coil2 : Coil
            The second coil.
        distance : float
            Axial distance between the coils in meters.

        Returns
        -------
        float
            The mutual inductance in Henries.

        Notes
        -----
        This is a lower-level method. For most use cases, prefer
        :meth:`calculate` which automatically selects the appropriate method.

        See Also
        --------
        calculate : Automatic method selection.
        calculate_near : Near-field calculation method.
        """
        coil_data = self._create_coil_data(coil1, coil2)
        precision_data = self._create_precision_data()

        return self.lib.lib.calculate_mutual_inductance_far(
            coil_data[0], precision_data[0], distance
        )

    def benchmark_near(self, n_repeats: int = 10000) -> None:
        """
        Run a benchmark for the near-field calculation method.

        Parameters
        ----------
        n_repeats : int, optional
            Number of iterations for the benchmark. Default is 10000.

        Notes
        -----
        Results are printed to stdout by the C library.
        """
        precision_data = self._create_precision_data()
        self.lib.lib.benchmark_mutual_inductance_near(precision_data[0], n_repeats)

    def benchmark_far(self, n_repeats: int = 10000) -> None:
        """
        Run a benchmark for the far-field calculation method.

        Parameters
        ----------
        n_repeats : int, optional
            Number of iterations for the benchmark. Default is 10000.

        Notes
        -----
        Results are printed to stdout by the C library.
        """
        precision_data = self._create_precision_data()
        self.lib.lib.benchmark_mutual_inductance_far(precision_data[0], n_repeats)

    def _create_coil_data(self, coil1: Coil, coil2: Coil) -> Any:
        """
        Create a CFFI CoilCalculationData struct.

        Parameters
        ----------
        coil1 : Coil
            The first coil.
        coil2 : Coil
            The second coil.

        Returns
        -------
        Any
            CFFI pointer to CoilCalculationData struct.
        """
        return self.lib.ffi.new(
            "CoilCalculationData *",
            {
                "r_1": coil1.inner_radius,
                "R_1": coil1.outer_radius,
                "L_1": coil1.length,
                "N_1": coil1.num_turns,
                "r_2": coil2.inner_radius,
                "R_2": coil2.outer_radius,
                "L_2": coil2.length,
                "N_2": coil2.num_turns,
            },
        )

    def _create_precision_data(self) -> Any:
        """
        Create a CFFI SumPrecisionData struct.

        Returns
        -------
        Any
            CFFI pointer to SumPrecisionData struct.
        """
        return self.lib.ffi.new(
            "SumPrecisionData *",
            {
                "k_terms": self.precision.k_terms,
                "l_terms": self.precision.l_terms,
                "n_terms": self.precision.n_terms,
            },
        )
