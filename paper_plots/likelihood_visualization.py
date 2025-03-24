import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl
import scipy as sc
from matplotlib.patches import ConnectionPatch
import matplotlib.ticker as ticker
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

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
facecol = 'lightgrey'

# Color definitions
rgba1 = (203 / 255, 139 / 255, 136 / 255, 255 / 255)
rgba2 = (136 / 255, 160 / 255, 203 / 255, 1)
rgba3 = (121 / 255, 192 / 255, 116 / 255, 1)

hbarc = 197.326
M = (938.272 + 939.565) / 2.0  # Averaged neutron/proton mass in MeV
units_factor = hbarc * hbarc / M

def Elab(p):
    return 2 * p ** 2 * hbarc ** 2 / M

def mom(E):
    return np.sqrt(M * E / 2 / hbarc ** 2)

# Load data
momentum_EMN_non_linear = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/momentum_EMN.txt')
energy_ = Elab(momentum_EMN_non_linear)
energy_EM = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/compare_phaseshifts/energy_EM.txt')
energy_ = energy_EM
phaseshifts_10010_old = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/phaseshift_files/phaseshifts_SVD/phaseshifts_unchanged_SLLJT_10010_lambda_1.80_s5.dat')

phaseshifts_10010 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/reference_phase_shifts/phaseshifts_unchanged_SLLJT_10010_lambda_1.80_s5.dat',skiprows=1)[:, 0]


#phaseshift_10010_random_sample = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/thesis_plots/files/phaseshifts_10010_1000000_to_100_3N_MCMC_25MeV.dat')
phaseshift_10010_random_sample = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/3N_MCMC/resampled_phaseshift_files/phaseshifts_10010_1000000_to_100_3N_MCMC_25MeV.dat')
phaseshift_10010_random_sample = phaseshift_10010_random_sample[1, :]

energy_lin = np.linspace(0.0007, 200, 2000)

energy_lin_likelihood = np.array([25, 50, 75, 100, 125, 150, 175, 200])

energy_EMN = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/energy_mesh_EMN_no_interpolation.txt')
#uncertainties_ = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/phaseshift_files/EKM_uncertainty/phaseshifts_uncertainties_SLLJT_10010_lambda_2.00.dat')[:, 3]
uncertainties_ = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/phaseshifts/EKM_uncertainties/phaseshifts_uncertainties_SLLJT_10010_lambda_2.00.dat', skiprows=1)[:, 3]


uncertainties_interpolate = sc.interpolate.interp1d(energy_EMN, uncertainties_)
uncertainties = uncertainties_interpolate(energy_lin)



ps_interpolate_old = sc.interpolate.interp1d(energy_, phaseshifts_10010_old)
phaseshifts_10010_old = ps_interpolate_old(energy_lin)


ps_interpolate = sc.interpolate.interp1d(energy_, phaseshifts_10010)
phaseshifts_10010 = ps_interpolate(energy_lin)

ps_interpolate_random_sample = sc.interpolate.interp1d(energy_, phaseshift_10010_random_sample)
phaseshifts_10010_random_sample = ps_interpolate_random_sample(energy_lin)


# Generate data for the main plot
x = energy_lin
y = phaseshifts_10010


linecol = 'dodgerblue'
errorcol = 'dimgrey'
onesigmacolor = 'chartreuse'




energy_lin_E1 = np.array([1, 3, 5, 10, 15, 20, 30, 50, 75, 100, 125, 150, 175, 200])
useless = np.ones((len(energy_lin_E1)))*50
# Create the main figure and plot
fig, ax = plt.subplots(figsize=(8, 6))
#ax.plot(energy_lin, phaseshifts_10010_old, color='red')



#ax.plot(x, phaseshifts_10010_random_sample, color='red', label='random sample')
ax.plot(x, y, label=r'$\delta_{\mathrm{ref}}$', linestyle='dotted', color='black')

#pwa93 = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/paper_plots/pwa93/3s1.txt')
#plt.plot(pwa93[:, 0], pwa93[:, 1])

ax.fill_between(x, y + uncertainties, y - uncertainties, color=facecol, alpha=1, label=r'N$^3$LO EKM')

ax.scatter(energy_lin_E1, useless, color='red', label = r'$\mathbf{E}_1$ likelihood grid', s=25)

for a, likelihood_energy in enumerate(energy_lin_likelihood):
    ax.axvline(likelihood_energy, ls='dashed', color=linecol ,label=r'$\mathbf{E}_2$ likelihood grid' if a == 0 else None)


#ax.yaxis.set_major_locator(ticker.MaxNLocator(10))
#ax.xaxis.set_major_locator(ticker.MaxNLocator(10))
ax.set_xticks(energy_lin_likelihood)
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
# Set labels and title for the main plot
#ax.set_title('Main plot with inset showing another dimension')
ax.set_xlabel(r'$E_{\mathrm{lab}}$ (MeV)')
ax.set_ylabel(r'$\delta$ (deg)')
ax.set_xlim(0, 200)
ax.set_ylim(0, 200)


# Mark the point at x=50 on the main plot
x_main = 50
y_main = ps_interpolate(x_main)
ax.plot(x_main, y_main, color=rgba2)  # Mark the point on the main plot


# Create the inset axes
# This places the inset in the top right corner
ax_inset = fig.add_axes([0.65, 0.45, 0.22, 0.35])

# Alternatively, use inset_axes for better control
# ax_inset = inset_axes(ax, width="30%", height="30%", loc='upper right')

# Plot data in the inset plot (another dimension)
x_norm = np.linspace(-3, 3, 100)
norm = sc.stats.norm(0, 1)
ax_inset.plot(x_norm, norm.pdf(x_norm), color='dimgrey')#, label='single energy likelihood')
ax_inset.axvline(0, color='black', ls='dotted')

ax_inset.axvspan(-1, 1, color=facecol)

# Remove the ticks on the inset plot
ax_inset.set_xticks([])
ax_inset.set_yticks([])
ax_inset.set_ylim(0, 0.42)
ax_inset.set_xlim(-3, 3)

# Optionally, set title or labels for the inset plot
ax_inset.set_title(r'$\mathcal{L}_{\delta(50\,\mathrm{MeV})}$', fontsize=18)

# Now, add the connection lines from the point at x=50 on the main plot to the inset plot

# Coordinates in the main plot
xy_main = (x_main, y_main)

# Coordinates in the inset plot (corners of the inset axes)
# Let's connect to the top-left and bottom-right corners of the inset plot
xy_inset_tl = (ax_inset.get_xlim()[0], ax_inset.get_ylim()[1])  # Top-left corner
xy_inset_br = (ax_inset.get_xlim()[1], ax_inset.get_ylim()[0])  # Bottom-right corner

# Create connection lines
con1 = ConnectionPatch(xyA=xy_main, coordsA=ax.transData,
                       xyB=xy_inset_tl, coordsB=ax_inset.transData,
                       axesA=ax, axesB=ax_inset, color="black", linestyle='--')
con2 = ConnectionPatch(xyA=xy_main, coordsA=ax.transData,
                       xyB=xy_inset_br, coordsB=ax_inset.transData,
                       axesA=ax, axesB=ax_inset, color="black", linestyle='--')

# Add the connection lines to the figure
fig.add_artist(con1)
fig.add_artist(con2)


#ax.legend(loc=2)
plt.figlegend(loc = (0.15, 0.55))
# Show the plot
plt.savefig('./plots/likelihood_visualization_2.pdf')
plt.show()
