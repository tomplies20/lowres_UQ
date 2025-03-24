import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib as mpl
import os
from phaseshift_calculator_LECs import *



output_folder = './phaseshift_percent_files'

singular_value_folder = '/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values'

#energy_ = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/energy_.txt')
energy_EM = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/compare_phaseshifts/energy_EM.txt')
energy_EMN = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/energy_mesh_EMN_no_interpolation.txt')

#EKM_uncertainty_folder = '/Users/pleazy/PycharmProjects/magic_quantification/library/phaseshift_files/EKM_uncertainty'
EKM_uncertainty_folder = '/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/phaseshift_files/EKM_uncertainty'
reference_folder = '/Users/pleazy/PycharmProjects/magic_quantification/library/phaseshift_files/phaseshifts_SVD'
#weights_folder = '/Users/pleazy/PycharmProjects/Phaseshift_Sampler/phaseshifts/weights'
#energy_ = np.loadtxt('/Users/pleazy/PycharmProjects/Phaseshift_Sampler/library/energy_.txt')

energy_lin_200 = np.linspace(energy_EM[0], 200, 200)
energy_lin = np.linspace(energy_EM[0], 200, 200)
lin_size = len(energy_lin)

partial_wave = '00001'
partial_wave = '10010'
#partial_wave = '01110'
#partial_wave = '11101'
#partial_wave = '11111'
#partial_wave = '11121'

SVD_rank = 4

EKM_uncertainty_file_name = f'phaseshifts_uncertainties_SLLJT_{partial_wave}_lambda_2.00.dat'
EKM_uncertainty_file = os.path.join(EKM_uncertainty_folder, EKM_uncertainty_file_name)
EKM_uncertainty = np.loadtxt(EKM_uncertainty_file)
N3LO_error = EKM_uncertainty[:, 3]
error_interpolate = sc.interpolate.interp1d(energy_EMN, N3LO_error)
error_interpolated = error_interpolate(energy_lin)

reference_file_name = f'phaseshifts_unchanged_SLLJT_{partial_wave}_lambda_1.80_s5.dat'
reference_file = os.path.join(reference_folder, reference_file_name)
reference_data = np.loadtxt(reference_file)
reference_N3LO = reference_data
reference_interpolate = sc.interpolate.interp1d(energy_EM, reference_N3LO)
reference_interpolated = reference_interpolate(energy_lin)





perc_steps = np.array([0.05, 0.1, 0.25, 0.5, 0.75, 1, -0.05, -0.1, -0.25, -0.5, -0.75, -1])
perc_steps = np.array([0.01, 0.025, 0.05, 0.1, 0.2, 0.5, -0.01, -0.025, -0.05, -0.1, -0.2, -0.5,])

#singular_value_file = f'SVD_chiral_order_N3LO_lambda_2.00_SLLJT_{partial_wave}_singular_values'
singular_value_file = f'singular_values_VNN_N3LO_EM500_SLLJT_{partial_wave}_lambda_1.80_Np_100_np_nocut.dat'
sv_path = os.path.join(singular_value_folder, singular_value_file)
# sv_path = f'/Users/pleazy/PycharmProjects/Phaseshift_Sampler/library/potentials/SVD_files/singular_values/SVD_chiral_order_N3LO_lambda_2.00_SLLJT_{partial_wave}_singular_values'
svs = np.loadtxt(sv_path)[0:SVD_rank+1]
phaseshifts_ = np.zeros([SVD_rank + 1, len(perc_steps), grid_size])
for i in range(SVD_rank+1):
    for j in range(len(perc_steps)):
        print(i)
        print(j)
        percent_plus = np.zeros([SVD_rank + 1])
        np.put(percent_plus, i, perc_steps[j])

        print(percent_plus)

        LECs = svs + svs*percent_plus
        phaseshifts_[i, j, :] = SVD(f'{partial_wave}', 3, SVD_rank, '1.80', LECs)

for u in range( SVD_rank+1):
    file_name = "phaseshifts_SLLJT_%s_lambda_1.80_s%s.dat" % (partial_wave, u + 1)
    output_file = os.path.join(output_folder, file_name)
    f = open(output_file, 'w')
    for m in range(grid_size):
        f.write(str(phaseshifts_[u, 0, m]) + ' ' + str(phaseshifts_[u, 1, m]) + ' ' + str(
            phaseshifts_[u, 2, m]) + ' ' + str(phaseshifts_[u, 3, m]) + ' ' + str(phaseshifts_[u, 4, m]) + ' ' +
                str(phaseshifts_[u, 5, m]) + ' ' + str(phaseshifts_[u, 6, m]) + ' ' + str(
            phaseshifts_[u, 7, m]) + ' ' + str(phaseshifts_[u, 8, m]) + ' ' + str(phaseshifts_[u, 9, m]) + ' ' +
                str(phaseshifts_[u, 10, m]) + ' ' + str(phaseshifts_[u, 11, m]) + "\n")






