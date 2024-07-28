import decimal
from decimal import Decimal, getcontext
import time

# Set the desired precision (e.g., 50 digits)
getcontext().prec = 30

max_n_terms = 100
max_l_terms = 100
max_k_terms = 100


class CoilCalculationData:
    def __init__(self, N_1, L_1, R_1, r_1, N_2, L_2, R_2, r_2):
        self.N_1 = Decimal(N_1)
        self.L_1 = Decimal(L_1)
        self.R_1 = Decimal(R_1)
        self.r_1 = Decimal(r_1)
        self.N_2 = Decimal(N_2)
        self.L_2 = Decimal(L_2)
        self.R_2 = Decimal(R_2)
        self.r_2 = Decimal(r_2)


def create_factorial_array(max_factorial):
    table = [Decimal(1)]
    for i in range(1, max_factorial + 1):
        table.append(table[-1] * Decimal(i))
    return table


def calculate_mutual_inductance_unoptimized(data, d, Z, timing=False):
    local_pi = Decimal('3.14159265358979323846264338327950288419716939937510')
    d = Decimal(d)
    Z = Decimal(Z)

    max_factorial = 2 * (max_l_terms - 1) + (max_n_terms - 1) + 2 * (max_k_terms - 1) + 1
    factorial_array = create_factorial_array(max_factorial)

    if timing:
        start_time = time.time()

    M_12 = Decimal(0)

    for n in range(0, max_n_terms):
        for l in range(0, max_l_terms):
            for k in range(0, max_k_terms):
                sign = Decimal(-1 if (l + k) % 2 == 1 else 1)

                numerator = (sign
                             * factorial_array[2 * l + 2 * k + n + 1]
                             * ((1 + data.L_1 / Z) ** (n + 1) - 1)
                             * ((data.R_1 / Z) ** (2 * l + 3) - (data.r_1 / Z) ** (2 * l + 3))
                             * ((data.R_2 / Z) ** (2 * k + 3) - (data.r_2 / Z) ** (2 * k + 3))
                             * (((1 + (data.L_1 + data.L_2 + d) / Z) ** (2 * l + 2 * k + n + 2)
                                 - (1 + (data.L_1 + d) / Z) ** (2 * l + 2 * k + n + 2)))
                             )
                denominator = (2 ** (2 * l + 2 * k + 2)
                               * Decimal((2 * l + 3) * (2 * k + 3))
                               * factorial_array[l] * factorial_array[l + 1]
                               * factorial_array[k] * factorial_array[k + 1]
                               * factorial_array[n + 1]
                               * ((1 + (data.L_1 + d) / Z) * (1 + (data.L_1 + data.L_2 + d) / Z)) ** (
                                       2 * l + 2 * k + n + 2)
                               )
                term = numerator / denominator

                if not term.is_finite():
                    break

                M_12 += term

    M_12 *= (Decimal(4e-7) * local_pi * local_pi * data.N_1 * data.N_2 * Z ** 5) \
            / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2))

    if timing:
        end_time = time.time()
        print('Elapsed time:', end_time - start_time)

    return M_12


def calculate_mutual_inductance(data, d, Z, timing=False):
    local_pi = Decimal('3.14159265358979323846264338327950288419716939937510')
    d_loc = Decimal(d)
    Z_loc = Decimal(Z)

    max_factorial = 2 * (max_l_terms - 1) + (max_n_terms - 1) + 2 * (max_k_terms - 1) + 1
    factorial_array = create_factorial_array(max_factorial)

    if timing:
        start_time = time.time()

    M_12 = Decimal(0)

    R_1_sq = data.R_1 * data.R_1
    r_1_sq = data.r_1 * data.r_1

    R_2_sq = data.R_2 * data.R_2
    r_2_sq = data.r_2 * data.r_2

    denom_1 = Decimal(1.0) / (Decimal(Z) + data.L_1 + Decimal(d))
    denom_2 = Decimal(1.0) / (Decimal(Z) + data.L_1 + data.L_2 + Decimal(d))

    denom_1_sq = denom_1 * denom_1
    denom_2_sq = denom_2 * denom_2

    L_1_plus_Z = data.L_1 + Z_loc

    loop_denom_1 = denom_1_sq
    loop_denom_2 = denom_2_sq
    power_2 = Decimal(4.0)

    loop_R_1 = R_1_sq * data.R_1
    loop_r_1 = r_1_sq * data.r_1

    for l in range(0, max_l_terms):
        if l > 0:
            loop_denom_1 *= denom_1_sq
            loop_denom_2 *= denom_2_sq
            power_2 *= Decimal(4.0)

            loop_R_1 *= R_1_sq
            loop_r_1 *= r_1_sq

        loop_R_1_sub_r_1 = loop_R_1 - loop_r_1

        save_first_loop_denom_1 = loop_denom_1
        save_first_loop_denom_2 = loop_denom_2
        save_first_power_2 = power_2

        loop_R_2 = R_2_sq * data.R_2
        loop_r_2 = r_2_sq * data.r_2

        for k in range(0, max_k_terms):

            if k > 0:
                loop_denom_1 *= denom_1_sq
                loop_denom_2 *= denom_2_sq
                power_2 *= Decimal(4.0)

                loop_R_2 *= R_2_sq
                loop_r_2 *= r_2_sq

            loop_R_2_sub_r_2 = loop_R_2 - loop_r_2
            sign = Decimal(-1 if (l + k) % 2 == 1 else 1)
            factor = Decimal((2 * l + 3) * (2 * k + 3))
            table_value = (factorial_array[l] * factorial_array[l + 1]
                           * factorial_array[k] * factorial_array[k + 1])
            loop_R_1_term_mul_R_2 = loop_R_1_sub_r_1 * loop_R_2_sub_r_2

            save_second_loop_denom_1 = loop_denom_1
            save_second_loop_denom_2 = loop_denom_2
            save_second_power_2 = power_2

            loop_L_1_plus_Z = L_1_plus_Z
            loop_Z = Z_loc

            for n in range(0, max_n_terms):

                if n > 0:
                    loop_denom_1 *= denom_1
                    loop_denom_2 *= denom_2

                    loop_L_1_plus_Z *= L_1_plus_Z
                    loop_Z *= Z_loc

                numerator = (sign
                             * factorial_array[2 * l + 2 * k + n + 1]
                             * loop_R_1_term_mul_R_2
                             * (loop_L_1_plus_Z - loop_Z)
                             * (loop_denom_1 - loop_denom_2)
                             )
                denominator = (power_2
                               * factor
                               * table_value
                               * factorial_array[n + 1]
                               )
                term = numerator / denominator

                if not term.is_finite():
                    break

                M_12 += term

            loop_denom_1 = save_second_loop_denom_1
            loop_denom_2 = save_second_loop_denom_2
            power_2 = save_second_power_2

        loop_denom_1 = save_first_loop_denom_1
        loop_denom_2 = save_first_loop_denom_2
        power_2 = save_first_power_2

    M_12 *= ((Decimal(4e-7) * local_pi * local_pi * Decimal(data.N_1 * data.N_2))
             / (data.L_1 * data.L_2 * (data.R_1 - data.r_1) * (data.R_2 - data.r_2))
             )

    if timing:
        end_time = time.time()
        print('Elapsed time:', end_time - start_time)

    return M_12


if __name__ == '__main__':
    data = CoilCalculationData(100, 0.1, 0.2, 0.1, 100, 0.1, 0.4, 0.3)
    d = 0.3
    Z = 1.2
    result = calculate_mutual_inductance_unoptimized(data, d, Z, timing=True)
    print(result)
