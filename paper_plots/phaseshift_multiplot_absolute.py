import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib as mpl
import os
import matplotlib.ticker as ticker

# ------------------------- Styling ------------------------------------------
softblack = 'k'  # Looks better when printed on tex file
gray = '0.7'

mpl.rcParams['figure.dpi'] = 180
mpl.rcParams['font.size'] = 25
mpl.rcParams['text.usetex'] = True
plt.rcParams["text.latex.preamble"] = r"\usepackage{lmodern}\usepackage{amsmath}"
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['axes.edgecolor'] = softblack
mpl.rcParams['axes.xmargin'] = 0
mpl.rcParams['axes.labelcolor'] = softblack
mpl.rcParams['axes.linewidth'] = 1.0

mpl.rcParams['axes.labelsize'] = 30  # Adjust the size as desired

# Set the tick label sizes
mpl.rcParams['xtick.labelsize'] = 30  # Adjust the size as desired
mpl.rcParams['ytick.labelsize'] = 30  # Adjust the size as desired



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
linecol='dodgerblue'
errorcol = 'dimgrey'
onesigmacolor = 'chartreuse'
#onesigmacolor = 'springgreen'


borders = np.loadtxt(
    '/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/phaseshift_distributions/nonplausible_borders.dat',
    skiprows=2
)

EKM_uncertainty_folder = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/EKM_uncertainties'
reference_folder        = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/reference_phase_shifts'

energy_EMN = np.loadtxt(
    '/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/energy_mesh_EMN_no_interpolation.txt'
)
energy_EM = np.loadtxt(
    '/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/compare_phaseshifts/energy_EM.txt'
)

energy_lin_200 = np.linspace(energy_EM[0], 200, 100)
energy_lin = np.linspace(energy_EM[0], 200, 100)
lin_size = len(energy_lin)

partial_waves = ['00001', '10010', '01110', '11101', '11111', '11121']
partial_wave_labels = [
    r'$^1$S$_0$', r'$^3$S$_1$', r'$^1$P$_1$',
    r'$^3$P$_0$', r'$^3$P$_1$', r'$^3$P$_2$'
]

data_folder = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/phaseshifts/resampled'



pwa_1s0 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/PWA93/1s0.txt')
pwa_3s1 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/PWA93/3s1.txt')
pwa_1p1 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/PWA93/1p1.txt')
pwa_3p0 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/PWA93/3p0.txt')
pwa_3p1 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/PWA93/3p1.txt')
pwa_3p2 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/PWA93/3p2.txt')
pwas = [pwa_1s0, pwa_3s1, pwa_1p1, pwa_3p0, pwa_3p1, pwa_3p2]
# ------------------------- Single-plot function ------------------------------
def single_plot(x, ax):
    """
    Creates a single subplot for partial_waves[x] on Axes ax.
    Leaves y tick labels on (no hiding).
    """
    partial_wave = partial_waves[x]
    partial_wave_label = partial_wave_labels[x]

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
        # 1) Interpolate from (energy_EM) to (energy_lin)
        ps_fun = sc.interpolate.interp1d(energy_EM, data[k, :])
        ps_temp = ps_fun(energy_lin)
        phaseshifts_interp[k] = ps_temp
        # 2) Re-interpolate to energy_lin_200, then back to energy_lin
        #ps_fun_2 = sc.interpolate.interp1d(energy_lin_200, ps_temp)
        #phaseshifts_interp[k] = ps_fun_2(energy_lin)

    # EKM band
    ax.fill_between(
        energy_lin,
        reference_interpolated - error_interpolated,
        reference_interpolated + error_interpolated,
        color=facecol, alpha=1.0, label=r'N$^3$LO EKM'
    )
    # ---------------- Plot all MCMC sample lines ----------------
    ax.plot(energy_lin, phaseshifts_interp[0],
            color=linecol, linewidth=1, alpha=0.5, label=r'$\mathbf{E}_2$ ppd')
    for irow in range(1, len(phaseshifts_interp)):
        ax.plot(energy_lin, phaseshifts_interp[irow],
                color=linecol, linewidth=1, alpha=0.5)

    # ---------------- 68% interval, EKM band, reference lines ----------------
    perc_16, perc_84 = np.percentile(phaseshifts_interp, [16, 84], axis=0)


    # 68% MCMC interval
    ax.plot(energy_lin, perc_84,  color=onesigmacolor, label=r'$68\%$ CI')
    ax.plot(energy_lin, perc_16,  color=onesigmacolor)



    # Reference line
    ax.plot(energy_lin, reference_interpolated, label=r'$\delta_{\mathrm{ref}}$',
            color='black', linestyle='dotted', linewidth=2)



    ax.plot(energy_lin, reference_interpolated + error_interpolated,
            color=errorcol, alpha=0.7)
    ax.plot(energy_lin, reference_interpolated - error_interpolated,
            color=errorcol, alpha=0.7)




    # ---------------- Labels & Title ----------------
    ax.text(
        0.6, 0.94, partial_wave_label, transform=ax.transAxes,
        fontsize=30, verticalalignment='top', horizontalalignment='left',
        bbox=dict(boxstyle='round', fc='white', ec='black', alpha=0.6)
    )

    #--------------------------
    #  plot pwa93 data
    #ax.plot(pwas[x][:, 0], pwas[x][:, 1], color='red')


    ax.set_xlabel(r'$E$ (MeV)')
    ax.set_ylabel(r'$\delta$ (deg)')
    ax.set_xlim(0, 200)

    # Tick locators
    ax.yaxis.set_major_locator(ticker.MaxNLocator(5))
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    # Turn on all side ticks
    ax.tick_params(axis="x", which="both", top=True, bottom=True)
    ax.tick_params(axis="y", which="both", left=True, right=True)

    # Hide y-label on middle and right columns
    if x not in [0, 3]:
        ax.set_ylabel('')

    # Hide x-label on top row
    if x not in [3, 4, 5]:
        ax.tick_params(axis="x", labelbottom=False)
        ax.set_xlabel('')


# ------------------------- Main multiplot function ---------------------------
def multiplot():
    fig, axes = plt.subplots(2, 3, figsize=(14, 10))

    # Create subplots for each partial wave
    for i in range(6):
        single_plot(i, axes[i // 3, i % 3])

    # Place the legend in the second subplot (axes[0,1])
    handles, labels = axes[0,1].get_legend_handles_labels()
    handles[0], handles[3] = handles[3], handles[0]
    labels[0], labels[3] = labels[3], labels[0]

    handles[0], handles[1] = handles[1], handles[0]
    labels[0], labels[1] = labels[1], labels[0]


    axes[0,1].legend(handles, labels, loc='center', bbox_to_anchor=(0.5, 0.54))
    # You can tweak loc or bbox_to_anchor for a different position

    # Adjust the spacing as needed
    fig.tight_layout(pad=1, w_pad=0.0, h_pad=1)

    # Save & show
    output_file = os.path.join(output_folder, 'phaseshifts_absolute_25MeV_resampled.pdf')
    plt.savefig(output_file)
    plt.show()


# ------------------------- Execute ------------------------------------------
if __name__ == "__main__":
    multiplot()
