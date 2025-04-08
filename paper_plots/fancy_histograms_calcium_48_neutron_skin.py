import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


# --------------------
# Setup Matplotlib RC
# --------------------
def setup_rc_params():
    softblack = 'k'  # Looks better when printed on tex file
    gray = '0.7'

    mpl.rcParams['figure.dpi'] = 180
    mpl.rcParams['font.size'] = 14
    mpl.rcParams['text.usetex'] = True
    plt.rcParams["text.latex.preamble"] = r"\usepackage{lmodern}\usepackage{amsmath}"
    mpl.rcParams['font.family'] = 'serif'
    mpl.rcParams['axes.edgecolor'] = softblack
    mpl.rcParams['axes.xmargin'] = 0
    mpl.rcParams['axes.labelcolor'] = softblack
    mpl.rcParams['axes.linewidth'] = 1.0
 # Adjust the size as desired

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
    #mpl.rcParams['xtick.major.width'] = 1.4
    #mpl.rcParams['ytick.major.width'] = 1.4
    # mpl.rcParams['xtick.major.width'] = 1
    # mpl.rcParams['ytick.major.width'] = 1

    mpl.rcParams['legend.edgecolor'] = 'inherit'  # inherits from axes.edgecolor
    mpl.rcParams['legend.facecolor'] = (1, 1, 1, 0.6)  # Set facecolor with alpha
    mpl.rcParams['legend.fancybox'] = True
    mpl.rcParams['legend.borderaxespad'] = 0.8
    mpl.rcParams['legend.framealpha'] = None  # facecolor alpha above
    mpl.rcParams['patch.linewidth'] = 0.8

    mpl.rc("savefig", transparent=False, bbox="tight", pad_inches=0.05, dpi=300, format="pdf")


# --------------------------------------------
# Function to plot histograms with summary stats
# --------------------------------------------
def plot_histograms(pdf_list, labels, colors, shifts=None, n_bins=10, baseline=0, separation=1, y_scale=1):
    """
    Plots vertically stacked histograms (using 10 common bins) with horizontal shifts.
    Also draws horizontal lines for the 68% (16th/84th) and 95% (2.5th/97.5th) percentiles,
    and marks the median with a white circle.

    Parameters:
      pdf_list: list of numpy arrays (each is a sample set, e.g. 99 data points)
      labels: list of labels for each histogram
      colors: list of colors for each histogram
      shifts: list of horizontal x-axis shifts for each histogram (if None, defaults to no shift)
      n_bins: number of common bins to use (default: 10)
      baseline: vertical starting value for the first histogram
      separation: vertical distance between histograms
    """
    if shifts is None:
        shifts = [0] * len(pdf_list)

    # Determine a common x-range from all data
    all_data = np.concatenate(pdf_list)
    xmin, xmax = np.min(all_data), np.max(all_data)
    # Create common bin edges and centers
    bin_edges = np.linspace(xmin, xmax, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    width = bin_edges[1] - bin_edges[0]

    fig, ax = plt.subplots(figsize=(5, 2 * len(pdf_list)))

    for i, data in enumerate(pdf_list):
        offset = baseline + i * separation
        shift = shifts[i]

        # Compute histogram (density normalized)
        counts, _ = np.histogram(data, bins=bin_edges, density=True)
        # Apply vertical scaling
        counts = counts * y_scale
        x_vals = bin_centers + shift
        ax.bar(x_vals, counts, width=width, bottom=offset, color=colors[i],
               edgecolor='gray', alpha=1, label=labels[i])

        # Compute summary statistics (percentiles)
        lower68 = np.percentile(data, 16) + shift
        upper68 = np.percentile(data, 84) + shift
        lower95 = np.percentile(data, 2.5) + shift
        upper95 = np.percentile(data, 97.5) + shift
        med = np.percentile(data, 50) + shift

        # Draw horizontal lines at the histogram baseline (offset)
        ax.hlines(offset, lower68, upper68, colors='gray', lw=6)
        ax.hlines(offset, lower95, upper95, colors='gray', lw=2)
        # Mark the median with a white circle
        ax.plot(med, offset, marker='o', color='white', markersize=6, zorder=10)

        print(f"{labels[i]}: median = {med - shift:.3f}, mean = {np.mean(data):.3f}, std = {np.std(data):.3f}")

    ax.set_xlabel(r"E($^{28}$O) (MeV)")
    ax.set_ylabel("Density + offset")
    ax.set_ylim(0, 2.4)
    ax.set_xlim(-185, -150)

    #ax.set_yticks([0.5, 1.5])
    #ax.set_xticklabels([r'E(})'])
    ax.set_yticks([])
    ax.set_ylabel("")
    ax.set_axisbelow(True)

    handles, labels = ax.get_legend_handles_labels()
    # Suppose you want to reverse the order:
    order = [1, 0]  # custom order indices for two handles
    ax.legend([handles[i] for i in order], [labels[i] for i in order], title=r'pr(E$(^{28}\mathrm{O})|\mathcal{D}$)')
    ax.grid(axis='x')
    plt.tight_layout()
    plt.savefig('./plots/histograms_O28_energies.pdf')
    plt.show()


# -----------------
# Example usage
# -----------------
setup_rc_params()

# Load your data (each file should have at least 99 rows in the second column)
pdf1 = np.loadtxt(
    '/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/energies_oxygen_28_25MeV_coupled.txt')[:, 1]
pdf2 = np.loadtxt(
    '/Users/pleazy/PycharmProjects/uncertainty_quantification/post_sampling/imsrg_data/energies_oxygen_28_1MeV_coupled.txt')[
       :, 1]
pdf_list = [pdf1, pdf2]
labels = [ r"$\mathbf{E}_2$", r"$\mathbf{E}_1$"]
colors = ["#aac7fc", "#f6b99b"]
shifts = [0, 0]  # Example: shift the second histogram by 0.2 along the x-axis

plot_histograms(pdf_list, labels, colors, shifts=shifts, n_bins=10, baseline=0.5, separation=1, y_scale=8)
