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

# Colors with adjusted alpha for transparency
rgba1 = (203 / 255, 139 / 255, 136 / 255, 0.7)
rgba2 = (136 / 255, 160 / 255, 203 / 255, 0.7)
rgba3 = (121 / 255, 192 / 255, 116 / 255, 0.7)

rgba1_line = (153 / 255, 85 / 255, 80 / 255, 1.0)  # Darker, less saturated red, fully opaque
rgba2_line = (60 / 255, 90 / 255, 160 / 255, 1.0)  # Darker, more saturated blue, fully opaque


# Load data for dataset 1
point_proton_radii_1 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/energies_oxygen_28_1MeV_coupled.txt')[:, 1]

neutron_skin_1 = point_proton_radii_1
percentile_1_16 = np.percentile(neutron_skin_1, 16)
percentile_1_84 = np.percentile(neutron_skin_1, 84)
percentile_1_50 = np.percentile(neutron_skin_1, 50)

# Load data for dataset 2
point_proton_radii_2 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/energies_oxygen_28_25MeV_coupled.txt')[:, 1]

print(np.mean(point_proton_radii_1), '+-', np.std(point_proton_radii_1))
print(np.mean(point_proton_radii_2), '+-', np.std(point_proton_radii_2))


neutron_skin_2 = point_proton_radii_2
percentile_2_16 = np.percentile(neutron_skin_2, 16)
percentile_2_84 = np.percentile(neutron_skin_2, 84)
percentile_2_50 = np.percentile(neutron_skin_2, 50)

# Define common bins based on both datasets
data_min = min(neutron_skin_1.min(), neutron_skin_2.min())
data_max = max(neutron_skin_1.max(), neutron_skin_2.max())
bins = np.linspace(data_min, data_max, 11)  # Adjust the number of bins as needed

# Create figure and axis with the same size as CODE_1
fig, ax = plt.subplots(figsize=(8, 5))



# Plot histograms
ax.hist(neutron_skin_1, bins=bins, color=rgba1, label=r'E$_1$ likelihood data', alpha=0.8, edgecolor='black')
ax.hist(neutron_skin_2, bins=bins, color=rgba2, label=r'E$_2$ likelihood data', alpha=0.7, edgecolor='black')


# Plot 68% confidence intervals as filled areas
ax.axvspan(percentile_1_16, percentile_1_84, edgecolor='r', facecolor=rgba1,  alpha=0.4, label=r'$68\%$ CI E$_1$ likelihood data',hatch='\\')

ax.axvspan(percentile_2_16, percentile_1_16, edgecolor='grey', facecolor=rgba2,  alpha=0.3, label=r'$68\%$ CI E$_2$ likelihood data',hatch='/')
ax.axvspan(percentile_1_84, percentile_2_84, edgecolor='grey', facecolor=rgba2,  alpha=0.3,hatch='/')

#ax.axvspan(percentile_1_16, percentile_1_84, color=rgba1,  alpha=0.4, label=r'68% CI E$_1$ likelihood data',hatch='-')
#ax.axvspan(percentile_2_16, percentile_2_84, color=rgba2,  alpha=0.3, label=r'68% CI E$_2$ likelihood data',hatch='|')




# Plot median lines
ax.axvline(percentile_1_50, color=rgba1_line, ls='-', lw=2)
ax.axvline(percentile_2_50, color=rgba2_line, ls='-', lw=2)

# Plot 68% lines
ax.axvline(percentile_1_84, color=rgba1_line, ls='--', lw=2)
ax.axvline(percentile_1_16, color=rgba1_line, ls='--', lw=2)

ax.axvline(percentile_2_84, color=rgba2_line, ls='--', lw=2)
ax.axvline(percentile_2_16, color=rgba2_line, ls='--', lw=2)


#non-implausible (nature PB208)
nonimp_color = 'orange'
#ax.axvline(-411.84 + 34.56, ls='dashed', color=nonimp_color, label='non-implausible 68% CI')
#ax.axvline(-411.84 - 34.56, ls='dashed', color=nonimp_color)
#ax.axvline(-411.84, ls='-', color=nonimp_color)
#ax.errorbar(x= -411.84, y=35, xerr=[[34.56], [34.56]], capsize=10, color=nonimp_color, fmt='o', markersize=10, label=r'non implausible $68\%$ CI')
#ax.set_ylim(0, 40)

ax.yaxis.set_visible(False)

# Remove duplicate labels in the legend
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
#ax.legend(by_label.values(), by_label.keys())

# Labels and title
ax.set_xlabel(r'$E(^{48}\mathrm{Ca})$ (MeV)')
#ax.set_xlim(-490, -370)
#ax.set_ylim(0, 50)
# ax.set_ylabel('Frequency')  # Uncomment if you wish to add a y-axis label

# Adjust layout and save the figure
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
#plt.tight_layout()
fig.subplots_adjust(left=0.05, right=0.95, bottom=0.175, top=0.95)
#fig.tight_layout(pad=1, w_pad=0.0, h_pad=1)
plt.savefig('./plots/energy_distribution_O28.pdf', dpi=300)
plt.show()
