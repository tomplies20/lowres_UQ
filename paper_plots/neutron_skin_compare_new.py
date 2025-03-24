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

#rgba2 = 'dodgerblue'

# Load data for dataset 1
point_proton_radii_1 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/Rp2c_calcium_48_1MeV_coupled.txt')[:, 1]
point_neutron_radii_1 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/R2nc_calcium_48_1MeV_coupled.txt')[:, 1]

neutron_skin_1 = np.sqrt(point_neutron_radii_1) - np.sqrt(point_proton_radii_1)
percentile_1_16 = np.percentile(neutron_skin_1, 16)
percentile_1_84 = np.percentile(neutron_skin_1, 84)
percentile_1_50 = np.percentile(neutron_skin_1, 50)

# Load data for dataset 2
point_proton_radii_2 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/Rp2c_calcium_48_25MeV_coupled.txt')[:, 1]
point_neutron_radii_2 = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/R2nc_calcium_48_25MeV_coupled.txt')[:, 1]

neutron_skin_2 = np.sqrt(point_neutron_radii_2) - np.sqrt(point_proton_radii_2)
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
ax.hist(neutron_skin_1, bins=bins, color=rgba1, label=r'$\mathbf{E}_1$ ppd', alpha=0.8, edgecolor='black')
ax.hist(neutron_skin_2, bins=bins, color=rgba2, label=r'$\mathbf{E}_2$ ppd', alpha=0.7, edgecolor='black')


# Plot 68% confidence intervals as filled areas
#ax.axvspan(percentile_1_16, percentile_1_84, color=rgba1, alpha=0.4, label=r'68% CI E$_1$ likelihood data', hatch='/')
#ax.axvspan(percentile_2_16, percentile_2_84, color=rgba2, alpha=0.3, label=r'68% CI E$_2$ likelihood data', hatch='/')

ax.axvspan(percentile_1_16, percentile_1_84, edgecolor='r', facecolor=rgba1,  alpha=0.4, label=r'$68\%$ CI $\mathbf{E}_1$ ppd',hatch='\\')

ax.axvspan(percentile_2_16, percentile_1_16, edgecolor='grey', facecolor=rgba2,  alpha=0.3, label=r'$68\%$ CI $\mathbf{E}_2$ ppd',hatch='/')
ax.axvspan(percentile_1_84, percentile_2_84, edgecolor='grey', facecolor=rgba2,  alpha=0.3,hatch='/')



# Plot median lines
ax.axvline(percentile_1_50, color=rgba1, ls='-', lw=2)
ax.axvline(percentile_2_50, color=rgba2, ls='-', lw=2)

#Plot 68% lines
ax.axvline(percentile_1_84, color=rgba1, ls='--', lw=2)
ax.axvline(percentile_1_16, color=rgba1, ls='--', lw=2)

ax.axvline(percentile_2_84, color=rgba2, ls='--', lw=2)
ax.axvline(percentile_2_16, color=rgba2, ls='--', lw=2)

ax.yaxis.set_visible(False)

#non-implausible (nature PB208)
nonimp_color = 'orange'
#ax.axvline(0.187, ls='dashed', color=nonimp_color, label='non-implausible 68% CI')
#ax.axvline(0.141, ls='dashed', color=nonimp_color)
#ax.axvline(0.164, ls='-', color=nonimp_color)
ax.errorbar(x= 0.164, y=35, xerr=[[0.164-0.141], [0.187-0.164]], capsize=10, color=nonimp_color, fmt='o', markersize=10, label=r'nonimplausible $68\%$ CI')
###gaute (nature 2016)
gaute_color='green'
#ax.axvline(0.1437, ls ='--', color=gaute_color, label='Hagen et al. magic results')
#ax.axvline(0.1426, ls ='--', color=gaute_color)
#ax.axvline(0.145, ls ='--', color=gaute_color)
#ax.axvline(0.1506, ls ='--', color=gaute_color)
#ax.axvline(0.14, ls ='--', color=gaute_color)

gaute_values = [0.1437, 0.1426, 0.145, 0.1506, 0.14]
percentile_gaute = np.percentile(gaute_values, [16, 50, 84])
gaute_mean = np.mean(gaute_values)
ax.errorbar(x= gaute_mean, y=32, xerr=[[gaute_mean -0.14], [0.1506 - gaute_mean]], capsize=10, color=gaute_color, fmt='x', markersize=10,label='Hagen et al. CCSD results')


ax.scatter(0.1294, 30, ls='--', color='red', label=r'Hagen et al. N$^2$LO$_{\mathrm{sat}}$ result', marker='v')

'''
# Remove duplicate labels in the legend
handles, labels = ax.get_legend_handles_labels()
by_label = dict(zip(labels, handles))
ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(0.42, 0.85))
#plt.figlegend(loc='center', ncol=1, labelspacing=1, bbox_to_anchor=(0.84, 0.28), fontsize=17)
'''
handles, labels = ax.get_legend_handles_labels()

# Define the order you want for the legend labels:
desired_order = [
    r'$\mathbf{E}_1$ ppd',
    r'$\mathbf{E}_2$ ppd',
    r'$68\%$ CI $\mathbf{E}_1$ ppd',
    r'$68\%$ CI $\mathbf{E}_2$ ppd',
    r'nonimplausible $68\%$ CI',
    'Hagen et al. CCSD results',
    r'Hagen et al. N$^2$LO$_{\mathrm{sat}}$ result'
]

# Build new lists of handles/labels in that order
ordered_handles = []
ordered_labels = []
for label in desired_order:
    if label in labels:
        idx = labels.index(label)
        ordered_handles.append(handles[idx])
        ordered_labels.append(labels[idx])

# Create the legend in the desired order

ax.legend(
    ordered_handles, ordered_labels,
    bbox_to_anchor=(0.395, 0.9),
    loc='best'
    #loc = 'upper left'
)




# Labels and title
ax.set_xlabel(r'$\Delta r_{\mathrm{np}}(^{48}\mathrm{Ca})$ (fm)')
# ax.set_ylabel('Frequency')  # Uncomment if you wish to add a y-axis label
ax.set_xlim(0.125, 0.19)
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# Adjust layout and save the figure
#plt.tight_layout()
#fig.tight_layout(pad=1, w_pad=0.0, h_pad=1)
#ax.set_position([0.15, 0.15, 0.8, 0.8])
fig.subplots_adjust(left=0.05, right=0.95, bottom=0.175, top=0.95)
plt.savefig('./plots/neutron_skin_distribution_ca48.pdf', dpi=300)
plt.show()
