from decimal import Decimal, getcontext
getcontext().prec = 25


def create_factorial_array(max_factorial):
    table = [Decimal(1)]
    for i in range(1, max_factorial + 1):
        table.append(table[-1] * Decimal(i))
    return table


def generate_lookup_table_near(num_l_terms, num_k_terms, num_n_terms):

    max_factorial = 2 * (num_l_terms - 1) + 2 * (num_k_terms - 1) + (num_n_terms - 1) + 1
    factorial_array = create_factorial_array(max_factorial)

    lookup_table = []
    for l in range(0, num_l_terms):
        local_matrix = []
        for k in range(0, num_k_terms):
            local_list = []
            for n in range(0, num_n_terms):
                sign = Decimal(-1 if (l + k) % 2 == 1 else 1)

                factor = Decimal((2 * l + 3) * (2 * k + 3))
                power_2 = Decimal(2.0) ** (2 * l + 2 * k + 2)
                leading_factorial = factorial_array[2 * l + 2 * k + n + 1]
                table_value = sign * leading_factorial / (
                        factor
                        * power_2
                        * factorial_array[l] * factorial_array[l + 1]
                        * factorial_array[k] * factorial_array[k + 1]
                        * factorial_array[n + 1])

                local_list.append(table_value)

            local_matrix.append(local_list)

        lookup_table.append(local_matrix)

    return lookup_table


def generate_lookup_table_far(num_l_terms, num_k_terms, num_n_terms):

    max_factorial = 2 * (num_l_terms - 1) + 2 * (num_k_terms - 1) + 2 * (num_n_terms - 1) + 1
    factorial_array = create_factorial_array(max_factorial)

    lookup_table = []
    for l in range(0, num_l_terms):
        local_matrix = []
        for k in range(0, num_k_terms):
            local_list = []
            for n in range(0, num_n_terms):
                sign = Decimal(-1 if (l + k) % 2 == 1 else 1)

                factor = Decimal((2 * l + 3) * (2 * k + 3))
                power_2 = Decimal(2.0) ** (2 * l + 2 * k + 2 * n + 2)
                leading_factorial = factorial_array[2 * l + 2 * k + 2 * n + 1]
                table_value = sign * leading_factorial / (
                        factor
                        * power_2
                        * factorial_array[l] * factorial_array[l + 1]
                        * factorial_array[k] * factorial_array[k + 1]
                        * factorial_array[2 * n + 1])

                local_list.append(table_value)

            local_matrix.append(local_list)

        lookup_table.append(local_matrix)

    return lookup_table


def write_lookup_table_to_file(lookup_table, max_terms, near: bool, fp_type="double"):
    filename = f"sum_lookup_table_{fp_type}_near.h" if near else f"sum_lookup_table_{fp_type}_far.h"
    table_name = "lookup_table_near" if near else "lookup_table_far"
    var_name = "MAX_TERMS_NEAR" if near else "MAX_TERMS_FAR"

    with open(filename, "w") as file:
        file.write("#include \"settings.h\"\n\n")
        file.write(f"#define {var_name} {max_terms}\n\n")
        file.write(f"static const {fp_type} {table_name}[{var_name}][{var_name}][{var_name}] =\n")
        file.write("{\n")
        for l in range(max_terms):
            file.write(f"  {{ // l = {l}\n")
            for k in range(max_terms):
                file.write(f"    {{ // k = {k}\n")  # Indent for clarity
                for n in range(max_terms):
                    if n == 0:
                        file.write("      ")
                    formatted_value = f"{lookup_table[l][k][n]:.20e}"
                    if (n + 1) % 4 == 0:
                        # End of a row of four, move to the next line
                        file.write(formatted_value + ",\n      ")
                    else:
                        file.write(formatted_value + ", ")
                file.write("},\n")
            file.write("  },\n")
        file.write("};\n")


if __name__ == '__main__':
    terms = 64
    table_near = generate_lookup_table_near(terms, terms, terms)
    write_lookup_table_to_file(table_near, terms, True, fp_type="double")

    table_far = generate_lookup_table_far(terms, terms, terms)
    write_lookup_table_to_file(table_far, terms, False, fp_type="double")
