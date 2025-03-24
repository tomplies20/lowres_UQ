import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker as ticker

'''
# Plot settings
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 1
plt.rcParams['ytick.major.width'] = 1
plt.rcParams['xtick.minor.size'] = 3
plt.rcParams['ytick.minor.size'] = 3
plt.rcParams['xtick.minor.width'] = 1
plt.rcParams['ytick.minor.width'] = 1
mpl.rcParams['axes.labelsize'] = 14  # Font size for x and y axis labels
mpl.rcParams['axes.titlesize'] = 12  # Font size for plot title
mpl.rcParams['legend.fontsize'] = 12  # Font size for legend
mpl.rcParams['xtick.labelsize'] = 12  # Font size for x-axis tick labels
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['font.family'] = 'DejaVu Sans'
'''

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


mpl.rcParams['axes.labelsize'] = 28  # Adjust the size as desired

# Set the tick label sizes
mpl.rcParams['xtick.labelsize'] = 28  # Adjust the size as desired
mpl.rcParams['ytick.labelsize'] = 28  # Adjust the size as desired



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
#mpl.rcParams['xtick.major.width'] = 1
#mpl.rcParams['ytick.major.width'] = 1

mpl.rcParams['legend.edgecolor'] = 'inherit'  # inherits from axes.edgecolor
mpl.rcParams['legend.facecolor'] = (1, 1, 1, 0.6)  # Set facecolor with alpha
mpl.rcParams['legend.fancybox'] = True
mpl.rcParams['legend.borderaxespad'] = 0.8
mpl.rcParams['legend.framealpha'] = None  # facecolor alpha above
mpl.rcParams['patch.linewidth'] = 0.8

text_bbox = dict(boxstyle='round', fc=(1, 1, 1, 0.6), ec=softblack, lw=0.8)

# Colors
rgba1 = (203 / 255, 139 / 255, 136 / 255, 0.7)  # Adjusted alpha for transparency
rgba2 = (136 / 255, 160 / 255, 203 / 255, 0.7)  # Adjusted alpha for transparency

rgba1_line = (153 / 255, 85 / 255, 80 / 255, 1.0)  # Darker, less saturated red, fully opaque
rgba2_line = (60 / 255, 90 / 255, 160 / 255, 1.0)  # Darker, more saturated blue, fully opaque

# Load data
#file_25MeV = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/oxygen/files/Ca48_Rp2c_1000000_to_100_samples_100_walkers_19_parameters_25MeV.txt')
Rp2c_25MeV = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/Rp2c_calcium_48_25MeV_coupled.txt')[:, 1]
Rso2_25MeV = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/Rso2_calcium_48_25MeV_coupled.txt')
radii_25MeV = np.sqrt(Rp2c_25MeV + 0.8409**2 - 28 / 20 * 0.1155 + 0.033 + Rso2_25MeV)

percentile_50_25MeV = np.percentile(radii_25MeV, 50)
percentile_16_25MeV = np.percentile(radii_25MeV, 16)
percentile_84_25MeV = np.percentile(radii_25MeV, 84)
#mean_25MeV = np.mean(radii_25MeV)

#file = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/oxygen/files/Ca48_Rp2c_1000000_to_100_samples_100_walkers_19_parameters_fixed_cov.txt')
Rp2c = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/Rp2c_calcium_48_1MeV_coupled.txt')[:, 1]
Rso2 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/Rso2_calcium_48_1MeV_coupled.txt')
radii = np.sqrt(Rp2c + 0.8409**2 - 28 / 20 * 0.1155 + 0.033 + Rso2)

print(np.std(radii_25MeV))

percentile_50 = np.percentile(radii, 50)
percentile_16 = np.percentile(radii, 16)
percentile_84 = np.percentile(radii, 84)
#mean = np.mean(radii)

# Define common bins
data_min = min(radii.min(), radii_25MeV.min())
data_max = max(radii.max(), radii_25MeV.max())
bins = np.linspace(data_min, data_max, 11)  # 20 bins

# Create figure and axis
fig, ax = plt.subplots(figsize=(8, 5))

# Plot histograms
ax.hist(radii, bins=bins, color=rgba1, label=r'E$_1$ likelihood data', alpha=0.8, edgecolor='black')
ax.hist(radii_25MeV, bins=bins, color=rgba2, label=r'E$_2$ likelihood data', alpha=0.7, edgecolor='black')

# Plot 68% confidence intervals with labels
#ax.axvspan(percentile_16, percentile_84, color=rgba1, alpha=0.3, label=r' 68% CI E$_1$ likelihood data', hatch='/')
#ax.axvspan(percentile_16_25MeV, percentile_84_25MeV, color=rgba2, alpha=0.3, label=r'68% CI E$_2$ likelihood data', hatch='/')

ax.axvspan(percentile_16, percentile_84, edgecolor='r', facecolor=rgba1,  alpha=0.4, label=r'68% CI E$_1$ likelihood data',hatch='\\')

ax.axvspan(percentile_16_25MeV, percentile_16, edgecolor='grey', facecolor=rgba2,  alpha=0.3, label=r'68% CI E$_2$ likelihood data',hatch='/')
ax.axvspan(percentile_84, percentile_84_25MeV, edgecolor='grey', facecolor=rgba2,  alpha=0.3,hatch='/')

# Plot mean lines with labels
ax.axvline(percentile_50, color=rgba1_line, ls='-', lw=2)#, label='Fixed Covariance Mean')
ax.axvline(percentile_50_25MeV, color=rgba2_line, ls='-', lw=2)#, label='25 MeV Mean')

#Plot 68% lines

ax.axvline(percentile_84, color=rgba1_line, ls='--', lw=2)
ax.axvline(percentile_16, color=rgba1_line, ls='--', lw=2)

ax.axvline(percentile_84_25MeV, color=rgba2_line, ls='--', lw=2)
ax.axvline(percentile_16_25MeV, color=rgba2_line, ls='--', lw=2)



#non-implausible (nature PB208)
nonimp_color = 'orange'
#ax.axvline(3.36 + 0.14, ls='--', color=nonimp_color, label='non-implausible 68% CI')
#ax.axvline(3.36 - 0.13, ls='--', color=nonimp_color)
#ax.axvline(3.36, ls='-', color=nonimp_color)

ax.errorbar(x= 3.36, y=45, xerr=[[0.13], [0.14]], capsize=10, color=nonimp_color, fmt='o', markersize=10, label='non-implausible 68% CI')


#ax.scatter(3.477, 40, marker='v',  label='Experimental value', color='black')

# Remove duplicate labels in the legend
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
#ax.legend(by_label.values(), by_label.keys())

# Labels and title
ax.set_xlabel(r'$R_{\mathrm{ch}}(^{48}\mathrm{Ca})$ (fm)')
#ax.set_ylabel('Frequency')
#ax.set_title('Histogram of Radii')

ax.yaxis.set_visible(False)
ax.set_xlim(3.05, 3.52)
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# Show plot
#plt.tight_layout()
#fig.tight_layout(pad=1, w_pad=0.0, h_pad=1)
fig.subplots_adjust(left=0.05, right=0.95, bottom=0.175, top=0.95)
plt.savefig('./plots/charge_radius_distribution_ca48.pdf', dpi=300)
plt.show()
