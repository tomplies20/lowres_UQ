import numpy as np
import os
from phaseshift_calculator_EMN import ps, ps_coupled




output_dir = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/EKM_uncertainties'




kmax = 6.0

# physical constants
hbarc = 197.326
M = (938.272 + 939.565) / 2.0
def momentum_scale(energy):
    pion_momentum = 139.57039 / hbarc * np.ones([len(energy)])

    Q = np.empty([len(energy)])
    for k in range(len(energy)):
        Q[k] = np.sqrt(M * energy[k] / (2 * hbarc * hbarc))

    lambda_b = 600
    Q = np.maximum.reduce([Q, pion_momentum]) / (lambda_b / hbarc)
    return Q


chiral_orders = ['LO', 'NLO', 'N2LO', 'N3LO', 'N4LO']

partial_waves = ['00001', '10010', '01110', '11101', '11111', '11121']




for partial_wave in partial_waves:
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


    output_filename = f'phaseshifts_uncertainties_SLLJT_{partial_wave}_lambda_2.00.dat'
    output_file = os.path.join(output_dir, output_filename)
    with open(output_file, 'w') as f:
        f.write('LO' + ' ' + 'NLO' + ' ' + 'N2LO' + 'N3LO' + 'N4LO' + 'Tlab')
        f.write('\n')
        for x in range(len(errors_order[0])):
            for y in range(5):
                f.write(str(errors_order[y,x]) + ' ' )
            f.write(str(energy[x]))
            f.write('\n')
