import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib as mpl
import os
import matplotlib.ticker as ticker

# ------------------------- Styling ------------------------------------------
softblack = 'k'
gray = '0.7'

mpl.rcParams['figure.dpi'] = 180
mpl.rcParams['font.size'] = 18
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
colors = [cmap(i / 5) for i in range(6)]
colors = [cmap(i / 5) for i in range(6)]
facecol = 'lightgrey'

# Paths & data
output_folder = './plots'

linecol = 'dodgerblue'
errorcol = 'dimgrey'
onesigmacolor = 'chartreuse'




EKM_uncertainty_folder = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/EKM_uncertainties'
reference_folder        = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/reference_phase_shifts'

energy_EMN = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/energy_mesh_EMN_no_interpolation.txt')
energy_EM  = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/compare_phaseshifts/energy_EM.txt')

energy_lin_200 = np.linspace(energy_EM[0], 200, 100)
energy_lin = np.linspace(energy_EM[0], 200, 100)
lin_size = len(energy_lin)

# --- Choose a single partial wave (example: ^1S0, which is 00001) ---
partial_wave = '00001'
partial_wave_label = r'$^1S_0$'

#partial_wave = '10010'
#partial_wave_label = r'$^3S_1$'

non_implausible_samples = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/nonimplausible/1s0.txt')
#non_implausible_samples = np.loadtxt('./non_implausible_samples_reordered/coupled_samples.txt')


data_folder = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/phaseshifts/resampled'


def plot_single_partial_wave(partial_wave, label):
    """
    Plots the absolute phaseshifts for one partial wave (e.g., ^1S0)
    using MCMC samples + reference + EKM uncertainty.
    """
    # ---------------- Load data for MCMC phaseshifts ----------------
    data_file_name = f'phaseshifts_{partial_wave}_1000000_to_100_3N_MCMC_25MeV_coupled.dat'
    data_file = os.path.join(data_folder, data_file_name)
    data = np.loadtxt(data_file)  # shape: (N, len(energy_EM))

    # ---------------- EKM uncertainty ----------------
    EKM_uncertainty_file_name = f'phaseshifts_uncertainties_SLLJT_{partial_wave}_lambda_2.00.dat'
    EKM_uncertainty_file = os.path.join(EKM_uncertainty_folder, EKM_uncertainty_file_name)
    EKM_uncertainty = np.loadtxt(EKM_uncertainty_file, skiprows=1)
    N3LO_error = EKM_uncertainty[:, 3]
    error_interpolate = sc.interpolate.interp1d(energy_EMN, N3LO_error)
    error_interpolated = error_interpolate(energy_lin)

    # ---------------- Reference data ----------------
    reference_file_name = f'phaseshifts_unchanged_SLLJT_{partial_wave}_lambda_1.80_s5.dat'
    reference_file = os.path.join(reference_folder, reference_file_name)
    reference_data = np.loadtxt(reference_file, skiprows=1)[:, 0]
    reference_interpolate = sc.interpolate.interp1d(energy_EM, reference_data)
    reference_interpolated = reference_interpolate(energy_lin)

    # ---------------- Interpolate MCMC samples to common energy grid ----------------
    phaseshifts_interp = np.zeros((len(data), lin_size))
    for k in range(len(data)):
        ps_fun = sc.interpolate.interp1d(energy_EM, data[k, :])
        ps_temp = ps_fun(energy_lin)
        phaseshifts_interp[k] = ps_temp

    # Create the figure & Axes
    fig, ax = plt.subplots(figsize=(8, 6))


    #Plot non implausible samples
    energy_grid_ni = [1, 5, 10, 25, 50, 100, 150, 200, 300]
    plt.plot(energy_grid_ni, non_implausible_samples[0, :], color='lightsalmon', alpha=0.6, linewidth=1, label='nonimplausible samples')
    for ps_1s0 in non_implausible_samples[1:]:
        plt.plot(energy_grid_ni, ps_1s0[:], color='lightsalmon', alpha = 0.6, linewidth=1)




    # Plot EKM band
    ax.fill_between(
        energy_lin,
        reference_interpolated - error_interpolated,
        reference_interpolated + error_interpolated,
        color=facecol, alpha=1.0
    )

    # Plot MCMC sample lines
    ax.plot(energy_lin, phaseshifts_interp[0],
            color=linecol, linewidth=1, alpha=0.5, label=r'$\mathbf{E}_2$ PPD')
    for irow in range(1, len(phaseshifts_interp)):
        ax.plot(energy_lin, phaseshifts_interp[irow],
                color=linecol, linewidth=1, alpha=0.5)

    # 68% MCMC interval
    perc_16, perc_84 = np.percentile(phaseshifts_interp, [16, 84], axis=0)
    ax.plot(energy_lin, perc_84, color=onesigmacolor, label=r'$68\%$ CI')
    ax.plot(energy_lin, perc_16, color=onesigmacolor)

    # Reference lines
    ax.plot(energy_lin, reference_interpolated,
            label=r'$\delta_{\mathrm{ref}}$', color='black',
            linestyle='dotted', linewidth=2)
    ax.plot(energy_lin, reference_interpolated + error_interpolated,
            color=errorcol, alpha=0.7)
    ax.plot(energy_lin, reference_interpolated - error_interpolated,
            color=errorcol, alpha=0.7, label=r'N$^3$LO EKM')


    '''
    pwa_ps = np.loadtxt('./pwa93/1s0.txt')
    energy = pwa_ps[:, 0]
    ps_1s0 = pwa_ps[:, 1]
    ax.plot(energy, ps_1s0, label='pwa93', color='red')
    '''


    # Title
    ax.text(0.85, 0.9, label, transform=ax.transAxes,
            fontsize=24, verticalalignment='top', horizontalalignment='left',
            bbox=dict(boxstyle='round', fc='white', ec='black', alpha=0.6))

    ax.set_xlabel(r'$E_{\mathrm{lab}}$ (MeV)')
    ax.set_ylabel(r'$\delta$ (deg)')
    ax.set_xlim(0, 200)
    ax.set_ylim(-25, 75)

    # Tick locators
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.tick_params(axis="x", which="both", top=True, bottom=True)
    ax.tick_params(axis="y", which="both", left=True, right=True)

    # Legend
    handles, labels = ax.get_legend_handles_labels()


    # Swap the first two elements.
    handles[0], handles[1] = handles[1], handles[0]
    labels[0], labels[1] = labels[1], labels[0]

    handles[3], handles[1] = handles[1], handles[3]
    labels[3], labels[1] = labels[1], labels[3]

    handles[3], handles[4] = handles[4], handles[3]
    labels[3], labels[4] = labels[4], labels[3]

    handles[2], handles[3] = handles[3], handles[2]
    labels[2], labels[3] = labels[3], labels[2]


    handles[3], handles[4] = handles[4], handles[3]
    labels[3], labels[4] = labels[4], labels[3]



    ax.legend(handles, labels, loc='lower left')

    # Save & show
    output_file = os.path.join(output_folder, f'phaseshift_{partial_wave}.pdf')
    #plt.tight_layout()
    plt.subplots_adjust(top=0.995, right=0.995, left=0.1, bottom=0.1)
    plt.savefig(output_file)
    plt.show()


# ------------------------- Execute the single plot --------------------------
if __name__ == "__main__":
    plot_single_partial_wave(partial_wave, partial_wave_label)
