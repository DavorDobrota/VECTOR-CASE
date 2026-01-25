from vector_case_py import (
    get_platform_info, Coil, SumPrecision, MutualInductanceCalculator
)


def main() -> None:
    """
    Demonstration of the VECTOR-CASE Python interface.

    This function shows basic usage of the library, including:

    - Creating coil objects
    - Calculating mutual inductance
    - Running benchmarks
    """
    print("=" * 60)
    print("VECTOR-CASE Python Interface Demo")
    print("=" * 60)
    print()

    # Show platform info
    info = get_platform_info()
    print("Platform Information:")
    for key, value in info.items():
        print(f"  {key}: {value}")
    print()

    # Define coil parameters
    coil1 = Coil(
        inner_radius=0.1,
        outer_radius=0.2,
        length=0.1,
        num_turns=100,
    )
    coil2 = Coil(
        inner_radius=0.3,
        outer_radius=0.4,
        length=0.1,
        num_turns=100,
    )
    coil3 = Coil(
        inner_radius=0.1,
        outer_radius=0.4,
        length=0.1,
        num_turns=100,
    )

    print(f"Coil 1: {coil1}\nCoil 2: {coil2}\nCoil 3: {coil3}\n")

    # Create calculator with default precision
    precision = SumPrecision(k_terms=32, l_terms=32, n_terms=32)
    calculator = MutualInductanceCalculator(precision=precision)

    print(f"Precision settings: {precision}\n")

    # Calculate mutual inductance
    distance = 0.1
    print(f"Mutual inductance at d = {distance} m for coils 1 and 3:")
    inductance = calculator.calculate(coil1, coil3, distance=distance, verbose=True)
    print(f"Mutual inductance (near):   {inductance:.10e} H\n")

    distance = 1.5
    print(f"Mutual inductance at d = {distance} m for coils 1 and 2:")
    inductance = calculator.calculate(coil1, coil2, distance=distance, verbose=True)
    print(f"Mutual inductance (far):    {inductance:.10e} H\n")

    distance = 0.1
    print(f"Mutual inductance at d = {distance} m for coils 1 and 2:")
    inductance = calculator.calculate(coil1, coil2, distance=distance, verbose=True)
    print(f"Mutual inductance (radial): {inductance:.10e} H\n")

    # Demonstrate near-field calculation with an explicit Z parameter
    print("For comparison, near field calculation with Z = 1.3:")
    z_param = 1.3
    distance = 0.1
    inductance_near = calculator.calculate_near(
        coil1, coil2, distance=distance, z_parameter=z_param
    )
    print(f"Mutual inductance (near):   {inductance_near:.10e} H\n")

    print(
        "Observe that the radial sum formula is more accurate than the "
        "near case, if applicable (outer radius of coil 2 is larger than "
        "inner radius of coil 1).\n"
    )

    # Run benchmarks
    calculator.benchmark_far(n_repeats=10000)
    print()

    calculator.benchmark_near(n_repeats=100)
    print()

    calculator.benchmark_radial(n_repeats=1000)
    print()

    calculator.benchmark_radial(n_repeats=100, fast=False)
    print()


if __name__ == "__main__":
    main()
