import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
import prettyplease as pp
from scipy.stats import norm, multivariate_normal
from matplotlib.patches import Ellipse
import matplotlib.ticker as ticker

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


rgba1 = (203 / 255, 139 / 255, 136 / 255, 255 / 255)

rgba2 = (136 / 255, 160 / 255, 203 / 255, 1)

rgba3 = (121 / 255, 192 / 255, 116 / 255, 1)

cmap = LinearSegmentedColormap.from_list(
    "custom_cmap",
    [(1, 1, 1, 1), rgba1]  # Light to custom color gradient
)

#samples = np.loadtxt('./resampled_parameters_MCMC_3N/3N_6000_no_cis_resampled_parameters.dat')
#samples = np.loadtxt('./resampled_parameters_MCMC_3N/3N_8000_1_resampled_parameters.dat')
#samples = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/phaseshift_distributions/files/MCMC_EKM_EM_3N_10000_no_cis_fixed_cov.txt')
samples = np.loadtxt('/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/samples/MCMC_EKM_EM_3N_10000_no_cis_25MeV_coupled_A.txt')
cD = samples[:, 17]
cE = samples[:, 18]

#cD = samples[:, 20]
#cE = samples[:, 21]

cDcE = np.array([cD, cE]).T
cov = np.cov(cDcE.T)
print(f'cov: {cov}')
mean = np.mean(cDcE.T, axis=1)
print(np.shape(cDcE))

col = linecol


coral = ['whitesmoke', rgba1]


cols = ['whitesmoke', col]

fig_x1x2, axes = pp.corner(cDcE, bins=30, labels=[r'$c_D$',r'$c_E$'], crosshairs=[1.264, -0.12],
                           quantiles=[0.16, 0.84], levels=(0.68,0.9),linewidth=1.0, plot_estimates=False,
                           colors=cols, n_uncertainty_digits=2,
                          title_loc='center',n_ticks=3 ,return_axes=True, figsize=(3,3))


fig_x1x2.subplots_adjust(
    left=0.22,   # more space on the left so "cD" isn't cut off
    bottom=0.21, # more space on the bottom so "cE" isn't cut off
    right=0.95, # adjust to taste
    top=0.87
    # optionally wspace=..., hspace=... for subplot spacing
)


ax = axes[1, 0]


'''
x_lim = ax.get_xlim()
y_lim = ax.get_ylim()

# Create a grid of points based on the x and y limits
x = np.linspace(x_lim[0], x_lim[1], 100)
y = np.linspace(y_lim[0], y_lim[1], 100)
x, y = np.meshgrid(x, y)
pos = np.dstack((x, y))
rv = multivariate_normal(mean, cov)
pdf = rv.pdf(pos)
ax.contour(x, y, pdf, levels=10, cmap='viridis')
'''

gauss_color = 'red'

data = cDcE
# Access the diagonal plots for x1 and x2
ax_x1 = axes[0, 0]  # The x1 histogram (top-left plot)
ax_x2 = axes[1, 1]  # The x2 histogram (bottom-right plot)


# Fit Gaussian for x1 and x2 (1D)
mean_x1, std_x1 = np.mean(data[:, 0]), np.std(data[:, 0])
mean_x2, std_x2 = np.mean(data[:, 1]), np.std(data[:, 1])

# Create a range of x values for plotting the Gaussian PDFs
x1_vals = np.linspace(*ax_x1.get_xlim(), 100)
x2_vals = np.linspace(*ax_x2.get_xlim(), 100)

# Extract the pre-existing histograms from the corner plot
hist_x1 = [child for child in ax_x1.get_children() if isinstance(child, plt.Polygon)][0]
hist_x2 = [child for child in ax_x2.get_children() if isinstance(child, plt.Polygon)][0]

# Get the vertices of the histogram to find the bin heights
verts_x1 = hist_x1.get_xy()
verts_x2 = hist_x2.get_xy()

# The y-values (heights) of the histogram bins are in the 2nd column of the vertices
bin_heights_x1 = verts_x1[:, 1]
bin_heights_x2 = verts_x2[:, 1]

# Find the maximum bin height for scaling the Gaussian to match the histogram
max_height_x1 = np.max(bin_heights_x1)
max_height_x2 = np.max(bin_heights_x2)

# Calculate the Gaussian PDF (unnormalized)
pdf_x1 = norm.pdf(x1_vals, mean_x1, std_x1)
pdf_x2 = norm.pdf(x2_vals, mean_x2, std_x2)

# Normalize the Gaussian PDFs to the peak of the histogram
pdf_x1 *= max_height_x1 / np.max(pdf_x1)
pdf_x2 *= max_height_x2 / np.max(pdf_x2)

# Plot the Gaussian fits on top of the PREEXISTING histograms (without re-plotting the histogram)
#ax_x1.plot(x1_vals, pdf_x1, color=gauss_color, lw=1.5)
#ax_x2.plot(x2_vals, pdf_x2, color=gauss_color, lw=1.5)#, label=r'1$\sigma$')

# Add legends to the diagonal plots
#ax_x1.legend()
#ax_x2.legend()
#plt.legend()

print(ax_x2.get_ylim())





# Access the 2D plot (bottom-left plot, i.e., `x1` vs. `x2`)
ax_2d = axes[1, 0]  # This accesses the (x1, x2) subplot in the corner plot

# Get the eigenvalues and eigenvectors of the covariance matrix
eigenvalues, eigenvectors = np.linalg.eigh(cov)

# Compute the angle of the ellipse from the largest eigenvector
angle = np.degrees(np.arctan2(*eigenvectors[:, 0][::-1]))

# Compute the 1-sigma and 2-sigma radii (scaling by the square root of the eigenvalues)
width_1sigma, height_1sigma = 2 * np.sqrt(eigenvalues)  # 1-sigma corresponds to sqrt of eigenvalues
width_2sigma, height_2sigma = 4 * np.sqrt(eigenvalues)  # 2-sigma corresponds to 2*sqrt of eigenvalues

# Add an ellipse for the 1-sigma region
ellipse_1sigma = Ellipse(mean, width_1sigma, height_1sigma, angle=angle,
                         edgecolor=gauss_color, fc='None', lw=1.5, label="1-Sigma")
#ax_2d.add_patch(ellipse_1sigma)

# Add an ellipse for the 2-sigma region
ellipse_2sigma = Ellipse(mean, width_2sigma, height_2sigma, angle=angle,
                         edgecolor=gauss_color, fc='None', lw=1.5, linestyle='--', label="2-Sigma")
# Add a scatter plot of the original data for reference
#ax_2d.scatter(data[:, 0], data[:, 1], s=2, c='blue', alpha=0.01)
#ax_2d.add_patch(ellipse_2sigma)



# Get the 16th and 84th percentiles for x1 and x2
x1_16th, x1_84th = np.percentile(data[:, 0], [16, 84])
x2_16th, x2_84th = np.percentile(data[:, 1], [16, 84])

# Plot vertical and horizontal lines indicating the 68% region
#ax_2d.axvline(x1_16th, color=rgba2, linestyle='--', label='68% Percentile (X1)')
#ax_2d.axvline(x1_84th, color=rgba2, linestyle='--')
#ax_2d.axhline(x2_16th, color=rgba2, linestyle='--', label='68% Percentile (X2)')
#ax_2d.axhline(x2_84th, color=rgba2, linestyle='--')


ax_x2.axvspan(x2_16th, x2_84th, hatch='/', alpha=0.3, color=col, label=r'$68\%$ CI')
ax_x1.axvspan(x1_16th, x1_84th, hatch='/', alpha=0.3, color=col)


# Determine the number of dimensions
ndim = data.shape[1]

# Loop over the axes to add minor ticks
for i in range(ndim):
    # Diagonal plots (histograms)
    ax = axes[i, i]
    # Get major ticks
    major_ticks = ax.get_xticks()
    # Compute minor ticks as midpoints between major ticks
    minor_ticks = (major_ticks[:-1] + major_ticks[1:]) / 2
    # Set minor ticks
    ax.set_xticks(minor_ticks, minor=True)
    # Enable minor ticks
    ax.minorticks_on()
    # Customize minor ticks appearance (optional)
   # ax.tick_params(axis='x', which='minor', length=3, width=1)

    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))

    # Lower triangle plots
    for j in range(i):
        ax = axes[i, j]
        # X-axis minor ticks
        major_ticks_x = ax.get_xticks()
        minor_ticks_x = (major_ticks_x[:-1] + major_ticks_x[1:]) / 2
        ax.set_xticks(minor_ticks_x, minor=True)
        ax.minorticks_on()
        #ax.tick_params(axis='x', which='minor', length=3, width=1)

        # Y-axis minor ticks
        major_ticks_y = ax.get_yticks()
        minor_ticks_y = (major_ticks_y[:-1] + major_ticks_y[1:]) / 2
        ax.set_yticks(minor_ticks_y, minor=True)
        ax.minorticks_on()
        #ax.tick_params(axis='y', which='minor', length=3, width=1)

        ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
        ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))



plt.savefig('./plots/cD_cE_corner_plot.pdf')
plt.show()



# Step 1: Extract the standard deviations (sqrt of diagonal elements)
std_devs = np.sqrt(np.diag(cov))

# Step 2: Calculate the correlation matrix
correlation_matrix = cov / np.outer(std_devs, std_devs)

# Step 3: Extract the off-diagonal element as the correlation coefficient
correlation_coefficient = correlation_matrix[0, 1]

print("Correlation Coefficient:", correlation_coefficient)

fig = plt.figure(figsize=(8, 6))
plt.scatter(cD, cE, color=rgba1)
plt.xlabel(r'$c_D$')
plt.ylabel(r'$c_E$')


#plt.savefig('cD_cE_corner_plot.pdf')


#plt.savefig('cD_cE_800000_to_100_samples_100_walkers_22_parameters.pdf')
#plt.show()