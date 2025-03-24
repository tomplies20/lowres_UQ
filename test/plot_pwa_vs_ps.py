import numpy as np
import matplotlib.pyplot as plt


pwa_3s1 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/PWA93/3p2.txt')

ps_3s1 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/reference_phase_shifts/phaseshifts_unchanged_SLLJT_11121_lambda_1.80_s5.dat', skiprows=1)


plt.plot(pwa_3s1[:, 0], pwa_3s1[:, 1])

plt.plot(ps_3s1[:, 1], ps_3s1[:, 0])


plt.xlim(0, 200)
plt.ylim(0, 20)
plt.show()