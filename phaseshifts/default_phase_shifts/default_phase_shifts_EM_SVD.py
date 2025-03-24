import numpy as np
import os
from phaseshift_calculator_LECs_EM import SVD, SVD_coupled



singular_value_dir = '/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values'

output_dir = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/reference_phase_shifts'

SVD_rank = 4

kmax = 6.5


partial_waves = ['00001', '10010', '01110', '11101', '11111', '11121']
#partial_waves = ['00001']
for partial_wave in partial_waves:
    singular_value_file_name = f'singular_values_VNN_N3LO_EM500_SLLJT_{partial_wave}_lambda_1.80_Np_100_np_nocut.dat'
    singular_value_file = os.path.join(singular_value_dir, singular_value_file_name)

    singular_values = np.loadtxt(singular_value_file)
    if partial_wave == '10010':
        phaseshifts, energy = SVD_coupled(partial_wave, '10210', '12010', '12210', SVD_rank, singular_values, kmax)
    if partial_wave == '11121':
        phaseshifts, energy = SVD_coupled(partial_wave, '11321', '13121', '13321', SVD_rank, singular_values, kmax)
    if partial_wave in ['00001', '01110', '11101', '11111']:
        print(partial_wave)
        phaseshifts, energy = SVD(partial_wave, SVD_rank,  singular_values, kmax)

    output_filename = f'phaseshifts_unchanged_SLLJT_{partial_wave}_lambda_1.80_s{str(SVD_rank+1)}.dat'
    output_file = os.path.join(output_dir, output_filename)
    with open(output_file, 'w') as f:
        f.write('N3LO_phase_shift Tlab')
        f.write('\n')
        for x, phaseshift in enumerate(phaseshifts):
            f.write(str(phaseshift) + ' ' + str(energy[x]))
            f.write('\n')