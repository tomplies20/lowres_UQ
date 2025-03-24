import os
import numpy as np
import scipy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator

# ------------------------- Global Styling -----------------------------------
softblack = 'k'  # Better for printing on TeX
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
mpl.rcParams['xtick.major.size'] = 7
mpl.rcParams['ytick.major.size'] = 7
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
colors = [cmap(i / 5) for i in range(6)]
facecol = 'lightgrey'

# ------------------------- Paths and Data -----------------------------------
output_folder = './plots'

# Example colors used
rgba1 = plt.get_cmap("copper")(0.8)
rgab1 = 'orange'
rgba2 = (136 / 255, 160 / 255, 203 / 255, 1)
rgba3 = (121 / 255, 192 / 255, 116 / 255, 1)

orange = '#E66100'
blue = '#5E9BD3'

rgba1 = 'indianred'
rgba2 = 'lightskyblue'
rgba2 = 'steelblue'

phaseshift_uncertainties_path = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/EKM_uncertainties'
phaseshift_reference_path = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/reference_phase_shifts'
data_folder = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/phaseshifts/full_set'
data_folder_resampled = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/phaseshifts/resampled'

# Load energies
energy_EMN = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/energy_mesh_EMN_no_interpolation.txt')
energy_EM = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/compare_phaseshifts/energy_EM.txt')

energy_lin_200 = np.linspace(energy_EM[0], 200, 200)
energy_lin = np.array([1, 5, 25, 50, 100, 150, 200])

partial_wave = '10010'

# Load MCMC data (25 MeV, resampled)
data_file_25 = os.path.join(data_folder_resampled,
                            f'phaseshifts_{partial_wave}_1000000_to_100_3N_MCMC_25MeV_coupled.dat')
mcmc_phaseshifts_all_25 = np.loadtxt(data_file_25)  # shape: (N, len(energy_EM))

# Load MCMC data (full)
data_file_full = os.path.join(data_folder,
                              f'phaseshifts_{partial_wave}_1000000_to_100_3N_MCMC_25MeV_coupled.dat')
mcmc_phaseshifts_all = np.loadtxt(data_file_full)

# Interpolation: 25 MeV (to 200 points, then to energy_lin)
mcmc_phaseshifts_interp_200_all_25 = np.zeros((len(mcmc_phaseshifts_all_25), 200))
for k in range(len(mcmc_phaseshifts_all_25)):
    f_25 = sc.interpolate.interp1d(energy_EM, mcmc_phaseshifts_all_25[k, :])
    mcmc_phaseshifts_interp_200_all_25[k] = f_25(energy_lin_200)

mcmc_phaseshifts_interp_all_25 = np.zeros((len(mcmc_phaseshifts_all_25), len(energy_lin)))
for k in range(len(mcmc_phaseshifts_all_25)):
    f_25b = sc.interpolate.interp1d(energy_lin_200, mcmc_phaseshifts_interp_200_all_25[k, :])
    mcmc_phaseshifts_interp_all_25[k] = f_25b(energy_lin)

# Interpolation: full MCMC (directly to energy_lin)
mcmc_phaseshifts_interp_all = np.zeros((len(mcmc_phaseshifts_all), len(energy_lin)))
for k in range(len(mcmc_phaseshifts_all)):
    f_full = sc.interpolate.interp1d(energy_EM, mcmc_phaseshifts_all[k, :])
    mcmc_phaseshifts_interp_all[k] = f_full(energy_lin)

# Load EKM uncertainties
EKM_file_name = f'phaseshifts_uncertainties_SLLJT_{partial_wave}_lambda_2.00.dat'
EKM_file = os.path.join(phaseshift_uncertainties_path, EKM_file_name)
EKM_data = np.loadtxt(EKM_file, skiprows=1)
N3LO_error = EKM_data[:, 3]

err_f = sc.interpolate.interp1d(energy_EMN, N3LO_error)
error_interp_200 = err_f(energy_lin_200)
err_f2 = sc.interpolate.interp1d(energy_lin_200, error_interp_200)
error_interpolated = err_f2(energy_lin)

# Load EKM reference phase shift
reference_file_name = f'phaseshifts_unchanged_SLLJT_{partial_wave}_lambda_1.80_s5.dat'
reference_file = os.path.join(phaseshift_reference_path, reference_file_name)
reference_data = np.loadtxt(reference_file, skiprows=1)[:, 0]

ref_f = sc.interpolate.interp1d(energy_EM, reference_data)
reference_200 = ref_f(energy_lin_200)
ref_f2 = sc.interpolate.interp1d(energy_lin_200, reference_200)
reference_interpolated = ref_f2(energy_lin)

# ------------------------- Plotting Functions -------------------------------
def single_plot(ax, i):
    """
    Plots a histogram for the i-th energy point, comparing
    the 25 MeV 'resampled' set and the full MCMC set.
    Also draws lines for EKM reference & 68% credible intervals.
    """
    data_25 = mcmc_phaseshifts_interp_all_25[:, i]
    data_full = mcmc_phaseshifts_interp_all[:, i]

    num_bins = 8
    bin_edges = np.linspace(
        min(np.min(data_25), np.min(data_full)),
        max(np.max(data_25), np.max(data_full)),
        num_bins + 1
    )

    # Energy label in top-right
    ax.text(0.965, 0.9, f'{energy_lin[i]} MeV',
            transform=ax.transAxes, ha='right', va='top',
            bbox=dict(boxstyle='round', fc='white', ec='black', alpha=0.6))

    # Plot histograms
    ax.hist(data_full, bins=bin_edges, alpha=0.6, edgecolor='dimgray',
            color=rgba1, density=True,
            label=r'$10^6$ ppd' if i == 6 else None)
    ax.hist(data_25, bins=bin_edges, alpha=0.5, edgecolor='dimgray',
            color=rgba2, density=True,
            label=r'$10^2$ ppd' if i == 6 else None)

    # Hide y-axis
    ax.axes.yaxis.set_visible(False)

    # Ticks
    ax.xaxis.set_major_locator(ticker.MaxNLocator(4))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(3))

    # EKM reference & uncertainty
    uncertainty = error_interpolated[i]
    ref_val = reference_interpolated[i]

    # 68% intervals
    pct_25 = np.percentile(data_25, [16, 50, 84])
    pct_full = np.percentile(data_full, [16, 50, 84])

    # Plot reference lines
    ax.axvline(ref_val, color='black',
               label=r'$\delta_{\mathrm{ref}}$' if i == 6 else None)
    ax.axvline(ref_val - uncertainty, color='black', linestyle='dashed',
               label=r'N$^3$LO EKM' if i == 6 else None)
    ax.axvline(ref_val + uncertainty, color='black', linestyle='dashed')

    # Plot 68% intervals
    ax.axvline(pct_full[0], color=rgba1, linestyle='dashed',
               label=r'$68\%$ CI $10^6$ ppd' if i == 6 else None)
    ax.axvline(pct_full[2], color=rgba1, linestyle='dashed')
    ax.axvline(pct_25[0], color=rgba2, linestyle='dashed',
               label=r'$68\%$ CI $10^2$ ppd' if i == 6 else None)
    ax.axvline(pct_25[2], color=rgba2, linestyle='dashed')

    # -- Larger x-limits logic (restored from your original approach) --
    # distance is the maximum of:
    #   (1) the EKM uncertainty (left & right),
    #   (2) the full 68% CI width from the large MCMC set,
    #   (3) the 68% CI width from the 25-MeV set.
    # Then we add extra padding (0.5 * distance * 2).
    # This yields broader x-limits as you had before.
    distance = 2 * np.max(np.array([
        uncertainty,
        pct_full[2] - pct_full[0],
        pct_25[2] - pct_25[0]
    ]))
    extra = 0.5 * distance * 2
    x_min = (ref_val - uncertainty) - extra
    x_max = (ref_val + uncertainty) + extra
    ax.set_xlim(x_min, x_max)

    return ax.get_legend_handles_labels()


def multiplot():
    """
    Creates a 2x4 grid of subplots. The last cell in the top row
    is used for the legend. The other 7 cells each show a histogram
    for one of the energy points in energy_lin.
    """
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))

    all_handles = []
    all_labels = []

    # We'll use 7 subplots for the 7 energies: (0,0),(0,1),(0,2),
    # (1,0),(1,1),(1,2),(1,3). The top-right cell (0,3) is purely for legend.
    idx_map = [(0, 0), (0, 1), (0, 2),
               (1, 0), (1, 1), (1, 2), (1, 3)]

    for subplot_idx, energy_idx in enumerate(range(len(energy_lin))):
        if subplot_idx == 7:
            break  # we only have 7 subplots for the data
        row, col = idx_map[subplot_idx]
        ax = axes[row, col]
        h, l = single_plot(ax, energy_idx)
        all_handles += h
        all_labels += l

    axes[1, 0].set_xlabel(r'$\delta$ (deg)')
    axes[1, 1].set_xlabel(r'$\delta$ (deg)')
    axes[1, 2].set_xlabel(r'$\delta$ (deg)')
    axes[1, 3].set_xlabel(r'$\delta$ (deg)')

    # Legend goes in top-right cell
    legend_ax = axes[0, 3]
    legend_ax.axis('off')

    # Deduplicate legend
    by_label = dict(zip(all_labels, all_handles))
    unique_handles = list(by_label.values())
    unique_labels = list(by_label.keys())
    legend_ax.legend(unique_handles, unique_labels, loc='center', frameon=False, fontsize=24)

    # Adjust spacing
    #fig.subplots_adjust(left=0.05, right=0.95, bottom=0.15, top=0.95, wspace=0.12) #wspace=0.08
    fig.subplots_adjust(left=0.005, right=0.995, bottom=0.1, top=0.995, wspace=0.12)
    out_file = os.path.join(output_folder,
                            f'phaseshift_distribution_25MeV_compare_{partial_wave}_subsampling.pdf')
    plt.savefig(out_file)
    plt.show()


if __name__ == "__main__":
    multiplot()
