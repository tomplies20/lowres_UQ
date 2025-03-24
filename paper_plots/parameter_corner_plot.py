import numpy as np
import emcee
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

import seaborn as sns
import prettyplease.prettyplease as pp
import matplotlib.ticker as ticker
import matplotlib as mpl

# ------------------------- Styling ------------------------------------------
softblack = 'k'  # Looks better when printed on tex file
gray = '0.7'


mpl.rcParams['figure.dpi'] = 180
mpl.rcParams['font.size'] = 20
mpl.rcParams['text.usetex'] = True
plt.rcParams["text.latex.preamble"] = r"\usepackage{lmodern}\usepackage{amsmath}"
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.edgecolor'] = softblack
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['axes.labelcolor'] = softblack
mpl.rcParams['axes.linewidth'] = 1.0

mpl.rcParams['ytick.direction'] = 'in'
mpl.rcParams['xtick.direction'] = 'in'
mpl.rcParams['xtick.color'] = softblack
mpl.rcParams['ytick.color'] = softblack
'''
mpl.rcParams['xtick.minor.size'] = 3.5
mpl.rcParams['ytick.minor.size'] = 3.5
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1.4
mpl.rcParams['ytick.major.width'] = 1.4
'''
mpl.rcParams['legend.edgecolor'] = 'inherit'
mpl.rcParams['legend.facecolor'] = (1, 1, 1, 0.6)
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.borderaxespad'] = 0.8
mpl.rcParams['legend.framealpha'] = None
mpl.rcParams['patch.linewidth'] = 0.8

linecol='dodgerblue'

###### import singular values


singular_value_ranges_percent_00001 = np.array([0.005, 0, 0.5, 0, 1]) # new 00001
singular_value_ranges_percent_10010 = np.array([0.02, 0.5, 1, 0, 0 ]) #10010
singular_value_ranges_percent_01110 = np.array([0.2, 0.05, 0.1, 0, 0]) # 01110
singular_value_ranges_percent_11101 = np.array([0.5, 0, 0, 0.1, 0]) #11101
singular_value_ranges_percent_11111 = np.array([0.05, 0, 0.05, 0.05, 0]) #11111
singular_value_ranges_percent_11121 = np.array([0.075, 0.15, 0.4, 0, 0]) #11121

svs_00001 = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values/singular_values_VNN_N3LO_EM500_SLLJT_00001_lambda_1.80_Np_100_np_nocut.dat')
svs_10010 = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values/singular_values_VNN_N3LO_EM500_SLLJT_10010_lambda_1.80_Np_100_np_nocut.dat')
svs_01110 = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values/singular_values_VNN_N3LO_EM500_SLLJT_01110_lambda_1.80_Np_100_np_nocut.dat')
svs_11101 = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values/singular_values_VNN_N3LO_EM500_SLLJT_11101_lambda_1.80_Np_100_np_nocut.dat')
svs_11111 = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values/singular_values_VNN_N3LO_EM500_SLLJT_11111_lambda_1.80_Np_100_np_nocut.dat')
svs_11121 = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values/singular_values_VNN_N3LO_EM500_SLLJT_11121_lambda_1.80_Np_100_np_nocut.dat')

all_singular_values = np.array([ svs_00001[0], svs_00001[2], svs_00001[4], svs_10010[0], svs_10010[1], svs_10010[2], svs_01110[0], svs_01110[1], svs_01110[2], svs_11101[0],svs_11101[1], svs_11101[3], svs_11111[0], svs_11111[2], svs_11111[3], svs_11121[0], svs_11121[1], svs_11121[2] ])



samples = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/samples/MCMC_EKM_EM_3N_10000_no_cis_25MeV_coupled_A.txt')
array = samples

# List to hold the separate arrays
separate_arrays = []

# Create 6 arrays, each starting from ith row and slicing every sixth row thereafter
num_walkers = 100
for i in range(num_walkers):
    sliced_array = array[i::num_walkers]  # Start from ith row, step by six
    separate_arrays.append(sliced_array)

print(np.shape(separate_arrays))



index_walker = 8

index_parameter = 1

param_1 = 15
param_2 = 16

#2, 10 trash
number_steps = 10000
number_parameters = 17 + 2

#burn_in = 2500
#burn_in = 1000
burn_in = 2000
burn_in = 0
#1, 3, 5
#1, 2, 3
#1, 2, 3
#1, 4
#1, 3, 4
#1, 2, 3
dict_LECs = {'c1_00001': 41.24475287917995, 'c2_00001': 3.6779529578686736, 'c3_00001': 2.240675329918265, 'c4_00001': 1.3154981170870734, 'c5_00001': 0.4105898445038659,
     'c1_10010': 56.73057943417189, 'c2_10010': 3.180500634833141, 'c3_10010': 0.7302127784460334, 'c4_10010': 0.3247620132567104, 'c5_10010': 0.2165281101925739,
             'c1_01110': 4.443979033264044, 'c2_01110': 2.9292881486682463, 'c3_01110': 1.469690501592981, 'c4_01110': 0.7123078990709594, 'c5_01110': 0.3375063600876752,
             'c1_11101': 3.026532911162979, 'c2_11101': 2.035025853323412, 'c3_11101': 0.5899304000517555, 'c4_11101': 0.4368381436558351, 'c5_11101': 0.128924046363905,
             'c1_11111': 3.968919734933947, 'c2_11111': 2.03651126236798, 'c3_11111': 1.6067962960276128, 'c4_11111': 0.6183216934167757, 'c5_11111': 0.2327552689700414,
             'c1_11121': 2.437107110171879, 'c2_11121': 1.001384771236426, 'c3_11121': 0.3870114142894703, 'c4_11121': 0.1450038365658212, 'c5_11121': 0.0516404886898715,
             'c1': -0.81, 'c3': -3.2, 'c4':  5.4, 'cD':1.264, 'cE':-0.12}
references = [41.24475287917995, 2.240675329918265, 0.4105898445038659,
              56.73057943417189, 3.180500634833141, 0.7302127784460334,
              4.443979033264044, 2.9292881486682463, 1.469690501592981,
              3.026532911162979, 0.4368381436558351,
              3.968919734933947, 1.6067962960276128, 0.6183216934167757,
              2.437107110171879, 1.001384771236426, 0.3870114142894703,
              1.264, -0.12]
labels = [r'$^{^1\mathrm{S}_0}s_1$', r'$^{^1\mathrm{S}_0}s_3$', r'$^{^1\mathrm{S}_0}s_5$',
          r'$^{^3\mathrm{S}_1}s_1$', r'$^{^3\mathrm{S}_1}s_2$', r'$^{^3\mathrm{S}_1}s_3$',
          r'$^{^1\mathrm{P}_1}s_1$', r'$^{^1\mathrm{P}_1}s_2$', r'$^{^1\mathrm{P}_1}s_3$',
          r'$^{^3\mathrm{P}_0}s_1$', r'$^{^3\mathrm{P}_0}s_4$',
          r'$^{^3\mathrm{P}_1}s_1$', r'$^{^3\mathrm{P}_1}s_3$', r'$^{^3\mathrm{P}_1}s_4$',
          r'$^{^3\mathrm{P}_2}s_1$', r'$^{^3\mathrm{P}_2}s_2$', r'$^{^3\mathrm{P}_2}s_3$',
          r'$c_D$', r'$c_E$']

separate_arrays = np.array(separate_arrays)
print(np.shape(separate_arrays))
separate_arrays = separate_arrays[:, burn_in:, :]
print(np.shape(separate_arrays))

one_walker = separate_arrays[index_walker]
#one_walker = one_walker[burn_in:]

print(np.shape(separate_arrays))
print(np.shape(one_walker))
#D1, D2, D3 = separate_arrays.shape
D1 = num_walkers
D2 = number_steps - burn_in
D3 = number_parameters
all_walkers = separate_arrays.reshape(D1*D2, D3)



singular_values_00001 = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values/singular_values_VNN_N3LO_EM500_SLLJT_00001_lambda_1.80_Np_100_np_nocut.dat')

colors = np.linspace(0, 1, len(one_walker))



# Sample data generation
np.random.seed(42)
k = 17 + 2 # Number of parameters
#data = np.random.randn(1000, k)  # 1000 samples of k parameters


col = linecol

cols = ['whitesmoke', col]

fig, axes = pp.corner(all_walkers, labels=labels, quantiles=[0.16, 0.84],crosshairs=references,  levels=(0.68,0.9),linewidth=1.0,
                          plot_estimates=False, colors=cols, n_uncertainty_digits=2,
                          title_loc='center', figsize=(8,8), return_axes=True, fontsize=13)

ndim=19
for col in range(ndim):
    ax_bottom = axes[-1, col]  # bottom row
    ax_left = axes[col, 0]
    ax_bottom.tick_params(axis='x', which='both', bottom=True, labelbottom=False)
    ax_left.tick_params(axis='y', which='both', left=True, labelleft=False)
    ax_left.yaxis.label.set_rotation(0)
    ax_left.yaxis.label.set_verticalalignment('center')
    ax_left.yaxis.label.set_horizontalalignment('right')

    ax_bottom.xaxis.label.set_rotation(45)


    # Grab the existing title, e.g. "cD\n1.23^{+0.12}_{-0.09}"
    old_title = axes[col, col].get_title()

    # If there is a newline, the second part after the split is the uncertainty
    splitted = old_title.split('\n', 1)
    if len(splitted) > 1:
        # splitted[0] is the variable name, splitted[1] is the "X +/- err"
        new_title = splitted[1]          # Keep only the second line
        axes[col, col].set_title(new_title,
                            x=2.15,    # move horizontally
                             y=1.4,   # move vertically
                             ha='right',
                             va='top', fontsize=11)  # Update the diagonal plot title
    else:                       #fontsize = 10 and x = 2.07 also works
        # No newline => nothing to remove
        pass
    #axes[col, col].set_title('')

    for a in range(ndim):
        for b in range (ndim):
            axes[a, b].tick_params(axis='x', which='both',
                       bottom=False, top=False,
                       labelbottom=False)
            axes[a, b].tick_params(axis='y', which='both',
                           left=False, right=False,
                           labelleft=False)

plt.savefig('./plots/corner_plot_E2.pdf')
plt.show()

