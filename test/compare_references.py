import numpy as np
import matplotlib.pyplot as plt

ps_ref_new = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/reference_phase_shifts/phaseshifts_unchanged_SLLJT_10010_lambda_1.80_s5.dat', skiprows=1)[:, 0]

ps_ref_old = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/phaseshift_files/phaseshifts_SVD/phaseshifts_unchanged_SLLJT_10010_lambda_1.80_s5.dat')

energy = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/reference_phase_shifts/phaseshifts_unchanged_SLLJT_10010_lambda_1.80_s5.dat', skiprows=1)[:, 1]


ps_pwa = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/PWA93/3s1.txt')

plt.plot(energy, ps_ref_old, label = 'old')
plt.plot(energy, ps_ref_new, label='new')
plt.plot(ps_pwa[:, 0], ps_pwa[:, 1])
plt.xlim(0, 100)
plt.legend()
plt.show()