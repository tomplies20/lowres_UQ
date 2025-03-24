import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import os
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as ticker
import matplotlib as mpl


# ------------------------- Styling ------------------------------------------
softblack = 'k'  # Looks better when printed on tex file
gray = '0.7'

mpl.rcParams['figure.dpi'] = 180
mpl.rcParams['font.size'] = 23
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
mpl.rcParams['xtick.minor.size'] = 3.5
mpl.rcParams['ytick.minor.size'] = 3.5
mpl.rcParams['xtick.minor.width'] = 1
mpl.rcParams['ytick.minor.width'] = 1
mpl.rcParams['xtick.major.size'] = 5
mpl.rcParams['ytick.major.size'] = 5
mpl.rcParams['xtick.major.width'] = 1.4
mpl.rcParams['ytick.major.width'] = 1.4

mpl.rcParams['legend.edgecolor'] = 'inherit'
mpl.rcParams['legend.facecolor'] = (1, 1, 1, 0.6)
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.borderaxespad'] = 0.8
mpl.rcParams['legend.framealpha'] = None
mpl.rcParams['patch.linewidth'] = 0.8

text_bbox = dict(boxstyle='round', fc=(1, 1, 1, 0.6), ec=softblack, lw=0.8)

# Example color setup
cmap = plt.get_cmap("copper")
colors = [cmap(i / 5) for i in range(6)]  # 6 steps (0 to 5)
facecol = 'lightgrey'

# ------------------------- Paths and Data -----------------------------------
output_folder = './plots'

rgba1 = (203 / 255, 139 / 255, 136 / 255, 255 / 255)
linecol = colors[4]
linecol = plt.get_cmap("copper")(0.8)
#linecol = 'lightsalmon'
#linecol='dodgerblue'
errorcol = 'dimgrey'
onesigmacolor = 'chartreuse'

rgba1 = (203 / 255, 139 / 255, 136 / 255, 255 / 255)

rgba1 = 'orange'
rgba1 = linecol

rgba2 = (136 / 255, 160 / 255, 203 / 255, 1)

rgba3 = (121 / 255, 192 / 255, 116 / 255, 1)



#LEC_path = '/Users/pleazy/PycharmProjects/Phaseshift_Sampler/phaseshifts/weighted_sampling_full/smarter_samples'
#phasehift_path = '/Users/pleazy/PycharmProjects/magic_quantification/3N_MCMC/resampled_phaseshift_files'
#####phaseshift_path = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/subsamples'


#phaseshifts_EM = '/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/phaseshift_distributions/files'
#####phaseshifts_EM = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/samples'

#phaseshifts_filtered_path = '/Users/pleazy/PycharmProjects/Phaseshift_Sampler/phaseshifts/weighted_sampling/data_filtered'

#phaseshifts_weighted_sampling_path = '/Users/pleazy/PycharmProjects/Phaseshift_Sampler/phaseshifts/weighted_sampling_full/smarter_samples'


#phaseshift_uncertainties_path = '/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/phaseshift_files/EKM_uncertainty'



#phaseshift_reference_path = '/Users/pleazy/PycharmProjects/magic_quantification/library/phaseshift_files/phaseshifts_SVD'


#weights_path = '/Users/pleazy/PycharmProjects/Phaseshift_Sampler/phaseshifts/initial_weights'



phaseshift_uncertainties_path = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/EKM_uncertainties'
phaseshift_reference_path = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/reference_phase_shifts'
energy_EMN = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/energy_mesh_EMN_no_interpolation.txt')
energy_EM = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/compare_phaseshifts/energy_EM.txt')

energy_lin_200 = np.linspace(energy_EM[0], 200, 200)
#energy_lin = np.linspace(energy_EM[0], 200, 200)
energy_lin = np.array([1, 3, 5, 10, 15, 20, 30, 50, 75, 100])
energy_lin = np.array([1, 5, 10, 20, 30, 50, 75, 100, 150, 200])
energy_lin = np.array([1, 5, 25, 50, 75, 100, 150, 200])
energy_lin = np.array([1, 5, 25, 50, 100, 150, 200])



#partial_wave = '00001'
partial_wave = '10010'
#partial_wave = '01110'
#partial_wave = '11101'
#partial_wave = '11111'
#partial_wave = '11121'


likelihood_25 = '25MeV'
likelihood_1 = 'fixed_cov'


lin_size = len(energy_lin)

dict_partial_wave = {
    '00001' : 0,
    '10010' : 1,
    '01110' : 2,
    '11101' : 3,
    '11111' : 4,
    '11121' : 5,

}


#100 resampled phase shifts, missleading variable name
data_folder = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/phaseshifts/full_set'
data_folder_resampled = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/phaseshifts/resampled'
# ---------------- Load data for MCMC phaseshifts ----------------
data_file_name_25 = f'phaseshifts_{partial_wave}_1000000_to_100_3N_MCMC_25MeV_coupled.dat'
data_file_25 = os.path.join(data_folder_resampled, data_file_name_25)
mcmc_phaseshifts_all_25 = np.loadtxt(data_file_25)  # shape: (N, len(energy_EM))


#data_file_name_1 = f'phaseshifts_{partial_wave}_1000000_to_100_3N_MCMC_25MeV_coupled.dat'
#data_file_1 = os.path.join(data_folder, data_file_name_1)
#mcmc_phaseshifts_all_1 = np.loadtxt(data_file_1)  # shape: (N, len(energy_EM))



# all 10^6 samples, no resampling

#mcmc_phaseshifts_all_1 = np.loadtxt(f'/Users/pleazy/PycharmProjects/magic_quantification/thesis_plots/files/phaseshifts_{partial_wave}_1000000_to_100_3N_MCMC_fixed_cov_X.dat')
#mcmc_phaseshifts_all_25 = np.loadtxt(f'/Users/pleazy/PycharmProjects/magic_quantification/thesis_plots/files/phaseshifts_{partial_wave}_1000000_to_100_3N_MCMC_25MeV.dat')








#mcmc_phaseshifts_all = np.loadtxt(f'/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/phaseshift_distributions/files/phaseshifts_EM_MCMC_{partial_wave}_3N_10000_no_cis_fixed_cov.dat')
#mcmc_phaseshifts_all = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/phaseshift_distributions/files/EKM_EM_phaseshift_3N_10000_no_cis.dat')[dict_partial_wave[partial_wave]::6, :]
#print(np.shape(mcmc_phaseshifts_all))
print(np.shape(energy_EM))
#phaseshifts_00001 = phaseshifts[::6, :]
#print(np.shape(mcmc_phaseshifts_all))
#plt.hist(mcmc_phaseshifts_all[:, 0], bins=20)
#plt.show()
mcmc_phaseshifts_interp_200_all_25 = np.zeros((len(mcmc_phaseshifts_all_25[:, 0]), 200))
for k in range(len(mcmc_phaseshifts_all_25[:, 0])):
    mcmc_ps_interpolate_200_all_25 = sc.interpolate.interp1d(energy_EM, mcmc_phaseshifts_all_25[k, :])
    mcmc_phaseshifts_interp_200_all_25[k] = mcmc_ps_interpolate_200_all_25(energy_lin_200)

mcmc_phaseshifts_interp_all_25 = np.zeros((len(mcmc_phaseshifts_all_25[:, 0]), lin_size))
for k in range(len(mcmc_phaseshifts_all_25[:, 0])):
    mcmc_ps_interpolate_all_25 = sc.interpolate.interp1d(energy_lin_200, mcmc_phaseshifts_interp_200_all_25[k, :])
    mcmc_phaseshifts_interp_all_25[k] = mcmc_ps_interpolate_all_25(energy_lin)


'''
mcmc_phaseshifts_interp_200_all_1 = np.zeros((len(mcmc_phaseshifts_all_1[:, 0]), 200))
for k in range(len(mcmc_phaseshifts_all_1[:, 0])):
    mcmc_ps_interpolate_200_all_1 = sc.interpolate.interp1d(energy_EM, mcmc_phaseshifts_all_1[k, :])
    mcmc_phaseshifts_interp_200_all_1[k] = mcmc_ps_interpolate_200_all_1(energy_lin_200)

mcmc_phaseshifts_interp_all_1 = np.zeros((len(mcmc_phaseshifts_all_1[:, 0]), lin_size))
for k in range(len(mcmc_phaseshifts_all_1[:, 0])):
    mcmc_ps_interpolate_all_1 = sc.interpolate.interp1d(energy_lin_200, mcmc_phaseshifts_interp_200_all_1[k, :])
    mcmc_phaseshifts_interp_all_1[k] = mcmc_ps_interpolate_all_1(energy_lin)
'''




#print(np.shape(mcmc_phaseshifts_all))









#mcmc_phaseshifts_1 = np.loadtxt(f'/Users/pleazy/PycharmProjects/magic_quantification/3N_MCMC/resampled_phaseshift_files/phaseshifts_{partial_wave}_1000000_to_100_3N_MCMC_fixed_cov.dat')#[:, :]
#mcmc_phaseshifts_25 = np.loadtxt(f'/Users/pleazy/PycharmProjects/magic_quantification/3N_MCMC/resampled_phaseshift_files/phaseshifts_{partial_wave}_1000000_to_100_3N_MCMC_25MeV.dat')#[:, :]





mcmc_phaseshifts_all = np.loadtxt(f'/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/phaseshifts/full_set/phaseshifts_10010_1000000_to_100_3N_MCMC_25MeV_coupled.dat')
#mcmc_phaseshifts_all = np.loadtxt(f'/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/phaseshift_distributions/files/phaseshifts_EM_MCMC_{partial_wave}_3N_10000_no_cis_fixed_cov.dat')
#mcmc_phaseshifts_all = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/phaseshift_distributions/files/EKM_EM_phaseshift_3N_10000_no_cis.dat')[dict_partial_wave[partial_wave]::6, :]
print(np.shape(mcmc_phaseshifts_all))
print(np.shape(energy_EM))
#phaseshifts_00001 = phaseshifts[::6, :]
print(np.shape(mcmc_phaseshifts_all))
#plt.hist(mcmc_phaseshifts_all[:, 0], bins=20)
#plt.show()
mcmc_phaseshifts_interp_all = np.zeros((len(mcmc_phaseshifts_all[:, 0]), len(energy_lin)))
for k in range(len(mcmc_phaseshifts_all[:, 0])):
    mcmc_ps_interpolate_200_all = sc.interpolate.interp1d(energy_EM, mcmc_phaseshifts_all[k, :])
    mcmc_phaseshifts_interp_all[k] = mcmc_ps_interpolate_200_all(energy_lin)





EKM_uncertainty_file_name = f'phaseshifts_uncertainties_SLLJT_{partial_wave}_lambda_2.00.dat'
EKM_uncertainty_file = os.path.join(phaseshift_uncertainties_path, EKM_uncertainty_file_name)
EKM_uncertainty = np.loadtxt(EKM_uncertainty_file, skiprows=1)
N3LO_error = EKM_uncertainty[:, 3]
error_interpolate = sc.interpolate.interp1d(energy_EMN, N3LO_error)
error_interpolated = error_interpolate(energy_lin_200)
error_interpolate_x = sc.interpolate.interp1d(energy_lin_200, error_interpolated)
error_interpolated = error_interpolate_x(energy_lin)

print('uncertainty')
print(error_interpolated)

reference_file_name = f'phaseshifts_unchanged_SLLJT_{partial_wave}_lambda_1.80_s5.dat'
reference_file = os.path.join(phaseshift_reference_path, reference_file_name)
reference_data = np.loadtxt(reference_file, skiprows=1)[:, 0]
reference_N3LO = reference_data

reference_interpolate = sc.interpolate.interp1d(energy_EM, reference_N3LO)
reference_interpolated = reference_interpolate(energy_lin_200)
reference_interpolate_x = sc.interpolate.interp1d(energy_lin_200, reference_interpolated)
reference_interpolated = reference_interpolate_x(energy_lin)
covariance_matrix = np.zeros((lin_size, lin_size))
signal_variance = np.abs(error_interpolated)
mean = reference_interpolated

f'phaseshifts_unchanged_SLLJT_{partial_wave}_lambda_2.00_s5.dat'






#titles = [r'$\delta$(x$\,$MeV)']
#tick_positions = [1,3, 5, 10, 15, 20, 30, 50, 75, 100]
tick_positions = energy_lin
tick_labels = [r'$\delta$({}$\,$MeV) (deg)'.format(x) for x in tick_positions]
#tick_positions = [r'$1$', r'$5$', r'$25$', r'$50$', r'$100$', r'$150$', r'$200$']
plot_labels = [r'${}\,$MeV'.format(x) for x in tick_positions]
the_label = r'$\delta$ (deg)'

def single_plot(ax, i):
    #ax = plt.subplot(2, 4, i+1)
    '''
    data1= phaseshifts_interp[:, i]
    data2 = phaseshifts_interp_filtered[:, i]
    data3 = phaseshifts_interp_weighted_sampling[:, i]
    '''
    #data_res = mcmc_phaseshifts_interp[:, i]
    data_res = data_25MeV = mcmc_phaseshifts_interp_all_25[:, i]
    #data_all = data_1MeV = mcmc_phaseshifts_interp_all_1[:, i]
    data_all = mcmc_phaseshifts_interp_all[:, i]
    #bin_edges = np.linspace(min(np.min(data1), np.min(data2), np.min(data3)),
     #                       max(np.max(data1), np.max(data2), np.max(data3)))

    num_bins = 8


    '''
    TEMPORARY FIX
    '''
    #data1 = data4

    # Calculate common bin edges based on the min and max of all datasets
    bin_edges = np.linspace(
        min(np.min(data_res), np.min(data_all)),
        np.max([np.max(data_res), np.max(data_all)]),
        num_bins + 1  # +1 because edges are one more than the number of bins
    )
    edges1 = bin_edges
    edges2 = bin_edges
    edges3 = bin_edges
    edges4 = bin_edges


    ax.text(
        0.9, 0.9, plot_labels[i], transform=ax.transAxes, verticalalignment='top', horizontalalignment='right',
        bbox=dict(boxstyle='round', fc='white', ec='black', alpha=0.6)
    )

    # Hide y-label on middle/right columns
    #if i not in [0, 3]:
    #    ax.set_ylabel('')
    # Hide x-label on top row
    #if i not in [3, 4, 5, 6]:
        #ax.tick_params(axis="x", labelbottom=False)
        #ax.set_xlabel('')

    # Plot the histograms
    if i == 6:
        #plt.hist(data2, bins=edges2, weights=weights_filtered, alpha=0.8, edgecolor='red', linestyle='dotted', histtype='step', linewidth=3, density=True, label = 'initial partial-wave distribution')
        #plt.hist(data3, bins=edges3, alpha=0.5, edgecolor='dimgray', color=rgba1,  density=True, label = '100 samples from partial-wave distribution')
        ax.hist(data_all, bins=bin_edges, alpha=0.6, edgecolor='dimgray', color=rgba1, density=True, label = r'$10^6$ ppd')
        ax.hist(data_res, bins = bin_edges, alpha=0.5, edgecolor='dimgray', color=rgba2, density=True, label = r'$10^2$ ppd')


    else:
        #plt.hist(data2, bins=edges2, weights=weights_filtered, alpha=0.8, edgecolor='red',
         #        linestyle='dotted', histtype='step', linewidth=3, density=True)
        #plt.hist(data3, bins=edges3, alpha=0.5, edgecolor='dimgray', color=rgba1, density=True)
        ax.hist(data_all, bins=bin_edges, alpha=0.6, edgecolor='dimgray', color=rgba1, density=True)
        ax.hist(data_res, bins=bin_edges, alpha=0.5, edgecolor='dimgray', color=rgba2, density=True)

    hist1, _ = np.histogram(data_all, bins=bin_edges, density=True)
    hist4, _ = np.histogram(data_res, bins=bin_edges, density=True)

    #plt.yticks([])
    max_1 = np.max(hist1)
    max_4 = np.max(hist4)
    max = np.max([max_1, max_4])
    print(max)
    old_y_ticks = np.linspace(0, max, 6)
    new_y_labels = [0,  .2,  .4,  .6,  .8,  1]
    plt.yticks(old_y_ticks, new_y_labels)

    #plt.gca().yaxis.set_minor_locator(MultipleLocator(0.2))
    #normalized_labels = [(x - 0) / (max - 0) for x in labels]

    ax.set_xlabel(tick_labels[i])
    ax.set_xlabel(the_label)

    #plt.gca().xaxis.set_major_locator(plt.MaxNLocator(4))

    ax.axes.yaxis.set_visible(False)
    ax.xaxis.set_major_locator(ticker.MaxNLocator(3))
    #ax.xaxis.set_major_locator(ticker.LinearLocator(3))
    #ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=1, integer=False, prune=None))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))





    uncertainty = error_interpolated[i]
    percentiles = np.percentile(data_res, [16, 50, 84])
    percentiles_all = np.percentile(data_all, [16, 50, 84])

    print(percentiles)
    print(percentiles_all)
    #percentiles_importance_resampling_samples = np.percentile(data3, [16, 50, 84])
    ####PREVIOUSLY data2 !!!!!!



    #plt.axvline(percentiles[1] + uncertainty, color='black')
    #plt.axvline(percentiles[1] - uncertainty, color='black')
    if i==6:
        ax.axvline(reference_interpolated[i] + uncertainty, color='black', linestyle='dashed', label = r' N$^3$LO EKM')
        ax.axvline(reference_interpolated[i] - uncertainty, color='black', linestyle='dashed')
        ax.axvline(reference_interpolated[i], color='black', label = r'$\delta_{\mathrm{ref}}$')


        ax.axvline(percentiles_all[0], color=rgba1, linestyle='dashed', label=r'$68\%$ CI $10^6$ ppd')
        ax.axvline(percentiles_all[2], color=rgba1, linestyle='dashed')

        ax.axvline(percentiles[0], color=rgba2, linestyle='dashed', label=r'$68\%$ CI $10^2$ ppd')
        ax.axvline(percentiles[2], color=rgba2, linestyle='dashed')


    else:
        ax.axvline(reference_interpolated[i] + uncertainty, color='black', linestyle='dashed')
        ax.axvline(reference_interpolated[i] - uncertainty, color='black', linestyle='dashed')
        ax.axvline(reference_interpolated[i], color='black')
        ax.axvline(percentiles[0], color=rgba2, linestyle='dashed')
        ax.axvline(percentiles[2], color=rgba2, linestyle='dashed')

        ax.axvline(percentiles_all[0], color=rgba1, linestyle='dashed')
        ax.axvline(percentiles_all[2], color=rgba1, linestyle='dashed')

    ##perhaps change this to weighted percentiles for the filtered distribution (or even the initial one)
    #plt.axvline(percentiles_importance_resampling_samples[0], color=rgba1, linestyle='dashed')
    #plt.axvline(percentiles_importance_resampling_samples[2], color=rgba1, linestyle='dashed')

    # if i % 4 == 0:
        #ax.ylabel('probability distribution')

    # Calculate the distance between the EKM uncertainty bounds

    distance = 2 * np.max(np.array([uncertainty,percentiles_all[2] - percentiles_all[0], percentiles[2] - percentiles[0]]))   # Since uncertainty is symmetric around reference_interpolated[i]

    # Calculate extra space (15% of the distance) to add on each side
    extra = 0.5 * distance *2

    # Set the new x-limits
    x_min = (reference_interpolated[i] - uncertainty) - extra
    x_max = (reference_interpolated[i] + uncertainty) + extra

    ax.set_xlim(x_min, x_max)

    if i in [0, 1, 2]:
        ax.set_xlabel('')

    return ax.get_legend_handles_labels()


def legend_plot(i):
    plt.subplot(2, 4, i)
    bbox_anchor = (0.7, 0.15)
    bbox_anchor = (0.87, 0.3)
    plt.figlegend(loc='center', ncol=1, labelspacing=1,  bbox_to_anchor=bbox_anchor, fontsize=20)
    plt.axis('off')

partial_waves = ['00001', '10010', '01110', '11101', '11111', '11121']



def multiplot():
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))

    # We'll plot in subplots 0..(lin_size-1), i.e. up to 6
    # The 8th cell (index=7) is reserved for legend
    #for idx in range(lin_size):


    all_handles = []
    all_labels  = []
    by_label = dict(zip(all_labels, all_handles))
    unique_handles = list(by_label.values())
    unique_labels  = list(by_label.keys())

    # Fill the first 5 cells with actual plots
    single_plot(axes[0, 0], 0)
    single_plot(axes[0, 1], 1)
    single_plot(axes[0, 2], 2)
    single_plot(axes[1, 0], 3)
    single_plot(axes[1, 1], 4)
    single_plot(axes[1, 2], 5)
    h, l, = single_plot(axes[1, 3], 6)

    all_handles.extend(h)
    all_labels.extend(l)


    # Legend subplot
    #legend_ax = ax_list[-1]
    legend_ax = axes[0, 3]
    # Hide the spines/ticks so we don't see a "box inside a box"
    #for spine in legend_ax.spines.values():
        #spine.set_visible(False)
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])

    # Deduplicate handles/labels
    by_label = dict(zip(all_labels, all_handles))
    unique_handles = list(by_label.values())
    unique_labels  = list(by_label.keys())

    # Place the legend in the center of that subplot, with no bounding box
    legend_ax.legend(unique_handles, unique_labels, loc='center', frameon=False)
    legend_ax.axis('off')
    #fig.tight_layout(pad=1, w_pad=0.0, h_pad=1)
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.15, top=0.95, wspace=0.08)
    #fig.tight_layout()
    out_file = os.path.join('./plots', f'phaseshift_distribution_25MeV_compare_{partial_wave}_subsampling.pdf')
    plt.savefig(out_file)
    plt.show()


if __name__ == "__main__":
    multiplot()


