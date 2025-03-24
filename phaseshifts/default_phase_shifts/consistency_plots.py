import numpy as np
import os
import matplotlib.pyplot as plt
from phaseshift_calculator_EMN import ps, ps_coupled
import scipy.interpolate as sc
from scipy.stats import norm



output_dir = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/EKM_uncertainties'




kmax = 6.0

# physical constants
hbarc = 197.326
M = (938.272 + 939.565) / 2.0
def momentum_scale(energy):
    pion_momentum = 139.57039 / hbarc * np.ones([len(energy)])
    #pion_momentum = 200 / hbarc * np.ones([len(energy)])

    Q = np.empty([len(energy)])
    for k in range(len(energy)):
        Q[k] = np.sqrt(M * energy[k] / (2 * hbarc * hbarc))

    lambda_b = 600
    Q = np.maximum.reduce([Q, pion_momentum]) / (lambda_b / hbarc)
    return Q


chiral_orders = ['LO', 'NLO', 'N2LO', 'N3LO', 'N4LO']

partial_waves = ['00001', '10010', '01110', '11101', '11111', '11121']


partial_waves = ['11101']


def compute(partial_wave):
    phaseshifts_chiral_orders = np.zeros((5, 100))
    for x, chiral_order in enumerate(chiral_orders):


        if partial_wave == '10010':
            phaseshifts, energy = ps_coupled(partial_wave, '10210', '12010', '12210', chiral_order, kmax)
        if partial_wave == '11121':
            phaseshifts, energy = ps_coupled(partial_wave, '11321', '13121', '13321', chiral_order,  kmax)
        if partial_wave in ['00001', '01110', '11101', '11111']:
            phaseshifts, energy = ps(partial_wave, chiral_order, kmax)
        phaseshifts_chiral_orders[x] = phaseshifts






    LO_phaseshifts = phaseshifts_chiral_orders[0]
    NLO_phaseshifts = phaseshifts_chiral_orders[1]
    N2LO_phaseshifts = phaseshifts_chiral_orders[2]
    N3LO_phaseshifts = phaseshifts_chiral_orders[3]
    N4LO_phaseshifts = phaseshifts_chiral_orders[4]


    diff_LO = np.abs(LO_phaseshifts)
    diff_NLO = np.abs(NLO_phaseshifts - LO_phaseshifts)
    diff_N2LO = np.abs(N2LO_phaseshifts - NLO_phaseshifts)
    diff_N3LO = np.abs(N3LO_phaseshifts - N2LO_phaseshifts)
    diff_N4LO = np.abs(N4LO_phaseshifts - N3LO_phaseshifts)


    Q = momentum_scale(energy)

    def max_LO():
        return Q * np.abs(LO_phaseshifts)


    def max_NLO():
        return np.maximum.reduce([Q ** 3 * diff_LO, Q ** 1 * diff_NLO])


    def max_N2LO():
        return np.maximum.reduce([Q ** 4 * diff_LO, Q ** 2 * diff_NLO, Q ** 1 * diff_N2LO])


    def max_N3LO():
        return np.maximum.reduce([Q ** 5 * diff_LO, Q ** 3 * diff_NLO, Q ** 2 * diff_N2LO, Q * diff_N3LO])


    def max_N4LO():
        return np.maximum.reduce(
            [Q ** 6 * diff_LO, Q ** 4 * diff_NLO, Q ** 3 * diff_N2LO, Q ** 2 * diff_N3LO, Q * diff_N4LO])


    # print(c_max_to_order(5))

    error_LO_x = max_LO()
    error_NLO_x = max_NLO()
    error_N2LO_x = max_N2LO()
    error_N3LO_x = max_N3LO()
    error_N4LO_x = max_N4LO()

    phaseshifts_order = [LO_phaseshifts, NLO_phaseshifts, N2LO_phaseshifts, N3LO_phaseshifts, N4LO_phaseshifts]

    errors_order = [error_LO_x, error_NLO_x, error_N2LO_x, error_N3LO_x, error_N3LO_x]
    errors_order = np.array(errors_order)
    print(np.shape(errors_order))


    expected_uncertainty_LO = Q * np.abs(LO_phaseshifts)
    expected_uncertainty_NLO = Q**3 * np.abs(LO_phaseshifts)
    expected_uncertainty_N2LO = Q**4 * np.abs(LO_phaseshifts)
    expected_uncertainty_N3LO = Q**5 * np.abs(LO_phaseshifts)
    expected_uncertainty_N4LO = Q**6 * np.abs(LO_phaseshifts)


    '''
    ratio_LO = diff_LO / expected_uncertainty_LO
    ratio_NLO = diff_NLO / expected_uncertainty_NLO
    ratio_N2LO = diff_N2LO / expected_uncertainty_N2LO
    ratio_N3LO = diff_N3LO / expected_uncertainty_N3LO
    ratio_N4LO = diff_N4LO / expected_uncertainty_N4LO
    '''

    ratio_LO = diff_LO / error_LO_x
    ratio_NLO = diff_NLO / error_NLO_x
    ratio_N2LO = diff_N2LO / error_N2LO_x
    ratio_N3LO = diff_N3LO / error_N3LO_x
    ratio_N4LO = diff_N4LO / error_N4LO_x


    num_points = 200

    energy_lin = np.linspace(energy[0], 200, num_points)


    interp_error_NLO = sc.interp1d(energy, error_NLO_x)
    error_NLO_interp = interp_error_NLO(energy_lin)

    interp_error_N2LO = sc.interp1d(energy, error_N2LO_x)
    error_N2LO_interp = interp_error_N2LO(energy_lin)

    interp_error_N3LO = sc.interp1d(energy, error_N3LO_x)
    error_N3LO_interp = interp_error_N3LO(energy_lin)

    interp_error_N4LO = sc.interp1d(energy, error_N4LO_x)
    error_N4LO_interp = interp_error_N4LO(energy_lin)


    interp_LO = sc.interp1d(energy, ratio_LO)
    ratio_LO_interp = interp_LO(energy_lin)

    interp_NLO = sc.interp1d(energy, ratio_NLO)
    ratio_NLO_interp = interp_NLO(energy_lin)

    interp_N2LO = sc.interp1d(energy, ratio_N2LO)
    ratio_N2LO_interp = interp_N2LO(energy_lin)

    interp_N3LO = sc.interp1d(energy, ratio_N3LO)
    ratio_N3LO_interp = interp_N3LO(energy_lin)

    interp_N4LO = sc.interp1d(energy, ratio_N4LO)
    ratio_N4LO_interp = interp_N4LO(energy_lin)

    # Define degrees of belief (probability thresholds)
    p_values = np.linspace(0.01, 0.99, 100)  # From 1% to 99%

    # Corresponding critical values from Gaussian distribution
    z_values = norm.ppf(0.5 + p_values / 2)  # Two-sided confidence intervals

    # Simulated observed differences (normally distributed for demonstration)
    #np.random.seed(42)
    #num_tests = 100
    #observed_diffs = np.random.normal(0, 1, num_tests)  # Simulated Delta X_n
    #estimated_errors = np.ones(num_tests)  # Assume δX_n = 1 for normalization

    # Compute success rate: fraction of cases within uncertainty bounds

    energy_bins = [(0, 50), (50, 100), (100, 200)]


    success_rates_NLO = []
    success_rates_N2LO = []
    success_rates_N3LO = []
    success_rates_N4LO = []

    for z in z_values:
        success_NLO = np.sum(np.abs(ratio_NLO_interp) <= z * error_NLO_interp) / 100
        success_rates_NLO.append(success_NLO)

        success_N2LO = np.sum(np.abs(ratio_N2LO_interp) <= z * error_N2LO_interp) / 100
        success_rates_N2LO.append(success_N2LO)

        success_N3LO = np.sum(np.abs(ratio_N3LO_interp) <= z * error_N3LO_interp) / 100
        success_rates_N3LO.append(success_N3LO)

        success_N4LO = np.sum(np.abs(ratio_N4LO_interp) <= z * error_N4LO_interp) / 100
        success_rates_N4LO.append(success_N4LO)

    colors = ['blue', 'green', 'red']
    #order = 'N2LO'
    for (E_min, E_max), color in zip(energy_bins, colors):
        # Select data in this energy range
        mask = (energy_lin >= E_min) & (energy_lin < E_max)
        selected_diffs = ratio_N3LO_interp[mask]
        selected_errors = error_N3LO_interp[mask]

        # Compute success rate for this energy interval
        success_rates = []
        #success_rate_std = []
        #success_rate_std_95 = []
        for z in z_values:
            success = np.sum(np.abs(selected_diffs) <= z * selected_errors) / len(selected_diffs)
            #success_std = np.sqrt(success * (1 - success) / len(selected_diffs))  # 68% band
            #success_std_95 = 2 * success_std  # 95% band (2σ)

            success_rates.append(success)
            #success_rate_std.append(success_std)
            #success_rate_std_95.append(success_std_95)

        success_rates = np.array(success_rates)
        #success_rate_std = np.array(success_rate_std)
        #success_rate_std_95 = np.array(success_rate_std_95)

        # Plot for this energy interval
        plt.plot(p_values, success_rates, label=f"{E_min}-{E_max} MeV", marker='o', color=color)

        # Add 95% confidence band (lighter shade)
        #plt.fill_between(p_values, success_rates - success_rate_std_95, success_rates + success_rate_std_95,
        #                 color=color, alpha=0.15, label=f"95% CI N3LO" if (E_min == 0) else "")

        # Add 68% confidence band (darker shade)
        #plt.fill_between(p_values, success_rates - success_rate_std, success_rates + success_rate_std,
        #                 color=color, alpha=0.3, label=f"68% CI N3LO" if (E_min == 0) else "")

    plt.text(0.7, 0.5, 'N3LO', fontsize=12, color='black',
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.4'))
    #plt.plot(p_values, success_rates_NLO, label='NLO')
    #plt.plot(p_values, success_rates_N2LO, label='N2LO')
    #plt.plot(p_values, success_rates_N3LO, label='N3LO')
    #plt.plot(p_values, success_rates_N4LO, label='N4LO')
    plt.plot(p_values, p_values, linestyle="dashed", color="black", label="Perfect Consistency (y=x)")
    plt.legend()
    plt.show()



for partial_wave in partial_waves:
    compute(partial_wave)


