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

    print(f"Coil 1: {coil1}")
    print(f"Coil 2: {coil2}")
    print()

    # Create calculator with default precision
    precision = SumPrecision(k_terms=32, l_terms=32, n_terms=64)
    calculator = MutualInductanceCalculator(precision=precision)

    print(f"Precision settings: {precision}")
    print()

    # Calculate mutual inductance
    distance = 0.1
    print(f"Calculating mutual inductance at distance d = {distance} m...")
    inductance = calculator.calculate(coil1, coil2, distance=distance, verbose=True)
    print(f"Mutual inductance: {inductance:.10e} H")
    print()

    # Demonstrate near-field calculation with an explicit Z parameter
    print("Near-field calculation with Z = 1.4:")
    z_param = 1.4
    inductance_near = calculator.calculate_near(
        coil1, coil2, distance=distance, z_parameter=z_param
    )
    print(f"Mutual inductance (near): {inductance_near:.10e} H")
    print()

    # Run benchmark
    print("Running near-field benchmark (10000 iterations)...")
    calculator.benchmark_near(n_repeats=10000)
    print()


if __name__ == "__main__":
    main()