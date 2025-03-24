import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import matplotlib as mpl
import os
from phaseshift_calculator_LECs import *





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
#plt.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 20  # Font size for x and y axis labels
mpl.rcParams['axes.titlesize'] = 14  # Font size for plot title
mpl.rcParams['legend.fontsize'] = 14  # Font size for legend
mpl.rcParams['xtick.labelsize'] = 16  # Font size for x-axis tick labels
mpl.rcParams['ytick.labelsize'] = 16
mpl.rcParams['font.family'] = 'DejaVu Sans'

rgba1 = (203 / 255, 139 / 255, 136 / 255, 255 / 255)

rgba2 = (136 / 255, 160 / 255, 203 / 255, 1)

rgba3 = (121 / 255, 192 / 255, 116 / 255, 1)

shades_of_grey = [
    (51/255, 51/255, 51/255),       # Dark Grey
    (102/255, 102/255, 102/255),    # Medium Dark Grey
    (153/255, 153/255, 153/255),    # Medium Grey
    (204/255, 204/255, 204/255),    # Medium Light Grey
    (230/255, 230/255, 230/255),    # Light Grey
    (242/255, 242/255, 242/255)     # Very Light Grey
]

colors = [
    (127/255, 179/255, 213/255),   # Pale Blue
    (157/255, 193/255, 131/255),   # Sage Green
    (177/255, 156/255, 217/255),   # Soft Lavender
    (216/255, 163/255, 180/255),   # Dusty Rose
    (233/255, 217/255, 133/255),   # Muted Yellow
    (204/255, 204/255, 204/255)    # Light Grey
]
#shades_of_grey = colors
shades_of_grey = shades_of_grey[::-1]

output_folder = './plots'

phaseshift_percent_folder = '/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/phaseshift_percent/phaseshift_percent_files'
EKM_uncertainty_folder = '/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/phaseshift_files/EKM_uncertainty'
reference_folder = '/Users/pleazy/PycharmProjects/magic_quantification/library/phaseshift_files/phaseshifts_SVD'

#energy_ = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/library/energy_.txt')
energy_EMN = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/uncertainties/energy_mesh_EMN_no_interpolation.txt')
energy_EM = np.loadtxt('/Users/pleazy/PycharmProjects/magic_quantification/phaseshifts/compare_phaseshifts/energy_EM.txt')
energy_lin_200 = np.linspace(energy_EM[0], 200, 200)
energy_lin = np.linspace(energy_EM[0], 200, 200)
lin_size = len(energy_lin)

partial_wave = '00001'
partial_wave = '10010'
partial_wave = '01110'
partial_wave = '11101'
#partial_wave = '11111'
#partial_wave = "11121"


dictionary_titles = {
    "00001": r'$^1$S$_0$',
    "10010": r'$^3$S$_1$',
    "01110": r'$^1$P$_1$',
    "11101": r'$^3$P$_0$',
    "11111": r'$^3$P$_1$',
    "11121": r'$^3$P$_2$'
}

x_max = 200

x_max_plot = 200

num_points = 200

SVD_rank = 7

def single_plot(partial_wave, i):
    plt.subplot(3, 3, i)

    file_name = "phaseshifts_SLLJT_%s_lambda_1.80_s%s.dat" % (partial_wave, i)
    phaseshift_file = os.path.join(phaseshift_percent_folder,file_name )

    phaseshifts_ = np.loadtxt(phaseshift_file)
    ps = np.zeros([num_points, len(phaseshifts_[1, :])])
    for k in range(len(phaseshifts_[1,:])):
        f_phaseshifts = sc.interpolate.interp1d(energy_EM, phaseshifts_[:,k])
        ps[:,k] = f_phaseshifts(energy_lin)



    EKM_uncertainty_file_name = f'phaseshifts_uncertainties_SLLJT_{partial_wave}_lambda_2.00.dat'
    EKM_uncertainty_file = os.path.join(EKM_uncertainty_folder, EKM_uncertainty_file_name)
    EKM_uncertainty = np.loadtxt(EKM_uncertainty_file)
    N3LO_error = EKM_uncertainty[:, 3]
    error_interpolate = sc.interpolate.interp1d(energy_EMN, N3LO_error)
    error_interpolated = error_interpolate(energy_lin)
    error_ = error_interpolated

    reference_file_name = f'phaseshifts_unchanged_SLLJT_{partial_wave}_lambda_1.80_s{SVD_rank +1}.dat'
    reference_file = os.path.join(reference_folder, reference_file_name)
    reference_data = np.loadtxt(reference_file)
    reference_N3LO = reference_data
    reference_interpolate = sc.interpolate.interp1d(energy_EM, reference_N3LO)
    reference_interpolated = reference_interpolate(energy_lin)
    reference = reference_interpolated

    grid_size = len(phaseshifts_[:, 11])




    #reference = phaseshifts[:, 10] #N3LO
    #error_ = phaseshifts[:, 11]

    alp = 1



    perc_steps = np.array([0.05, 0.1, 0.25, 0.5, 0.75, 1, -0.05, -0.1, -0.25, -0.5, -0.75, -1])
    perc_steps = np.array([0.01, 0.025, 0.05, 0.1, 0.2, 0.5, -0.01, -0.025, -0.05, -0.1, -0.2, -0.5, ])
    #perc_steps = np.array([0.5, 0.75, 1, 2, 5, 10, -0.5, -0.75, -1, -2, -5, -10])

    color_arr = ['grey', 'blue', 'purple', 'red', 'orange', 'yellow']
    color_arr = shades_of_grey
    if i==5:


        plt.plot(energy_lin, reference/reference - 1, label=r'$\delta_{\mathrm{ref}}$', color='black', linestyle='dotted', linewidth=2)
        for n in range(6):

            #plt.plot(energy_lin, ps[:, n])

            plt.plot(energy_lin,  ps[:,n]/reference - 1, label=str(perc_steps[n]) , color=color_arr[n], linewidth=2)  # 1 = NLO


        plt.fill_between(energy_lin, (reference + error_) / reference - 1,  (reference - error_) / reference - 1,
                         color=rgba1, alpha=alp,
                         label='error')

        for m in range(6):
            plt.plot(energy_lin, ps[:,m+6]/reference - 1, label=str(perc_steps[m+6]),color=color_arr[m], linewidth=2, linestyle='dashed')



    else:
        plt.plot(energy_lin,  reference / reference - 1, color='black',
                 linestyle='dotted', linewidth=2)
        plt.fill_between(energy_lin,  (reference + error_) / reference - 1,  (reference - error_) / reference - 1,
                         color=rgba1, alpha=alp,
                        )
        for n in range(6):
            # plt.plot(energy_lin, ps[:, n])

            plt.plot(energy_lin, ps[:, n] / reference - 1, color=color_arr[n],
                     linewidth=2)  # 1 = NLO

        for m in range(6):
            plt.plot(energy_lin,ps[:, m + 6] / reference - 1, color=color_arr[m],

                    linewidth=2, linestyle='dashed')



    plt.ylim(-0.1, 0.1)

    subplot_title = fr's$_{i}$'
    plt.text(0.05, 0.95, subplot_title, transform=plt.gca().transAxes, fontsize=22, verticalalignment='top', horizontalalignment='left', backgroundcolor='white')
    plt.xlabel('E (MeV)')
    plt.ylabel(r'$\Delta \delta_{\mathrm{rel}}$')
    plt.xlim(0, x_max_plot)

    #plt.legend(loc='center right')


def legend_plot(i):
    plt.subplot(3, 3, i)
    bbox_anchor = (0.84, 0.28)
    plt.figlegend(loc='center', ncol=2, labelspacing=1,  bbox_to_anchor=bbox_anchor, fontsize=20)
    plt.axis('off')
    handles, labels = plt.gca().get_legend_handles_labels()
    print(handles)
    print(labels)
    #plt.figlegend(handles, labels, loc='lower right', ncol=3, labelspacing=0.)


def multiplot(partial_wave):
    #plt.figure(figsize=(16, 8))
    plt.figure(figsize=(18, 18))

    single_plot(partial_wave, 1)
    single_plot(partial_wave, 2)
    single_plot(partial_wave, 3)
    single_plot(partial_wave, 4)
    single_plot(partial_wave, 5)
    single_plot(partial_wave, 6)
    single_plot(partial_wave, 7)
    single_plot(partial_wave, 8)
    legend_plot(9)

    #single_plot(partial_wave, 6)
    #single_plot(partial_wave, 7)
    #single_plot(partial_wave, 8)



    #plt.title(dictionary_titles[partial_wave])
    #plt.text(1, 1, dictionary_titles[partial_wave])
    plt.tight_layout()
    output_file = os.path.join(output_folder, partial_wave + '_percent_.pdf')

    #plt.savefig(output_file)
    plt.show()


multiplot(partial_wave)