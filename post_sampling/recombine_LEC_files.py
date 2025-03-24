import numpy as np
import os


output_path = 'singular_values_post'

resampled_LECs_path = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/subsamples'
#resampled_LECs_filename = '3N_10000_no_cis_resampled_parameters.dat'
resampled_LECs_filename = '3N_10000_no_cis_1MeV_coupled.dat'
resampled_LECs_file = os.path.join(resampled_LECs_path, resampled_LECs_filename)

resampled_LECs = np.loadtxt(resampled_LECs_file)


singular_values_path = '/Users/pleazy/PycharmProjects/magic_quantification/library/potentials/SVD_files_no_interpolation/singular_values'

partial_waves = ['00001', '10010', '01110', '11101', '11111', '11121']

singular_values_00001 = np.loadtxt(os.path.join(singular_values_path, f'singular_values_VNN_N3LO_EM500_SLLJT_00001_lambda_1.80_Np_100_np_nocut.dat'))
singular_values_10010 = np.loadtxt(os.path.join(singular_values_path, f'singular_values_VNN_N3LO_EM500_SLLJT_10010_lambda_1.80_Np_100_np_nocut.dat'))
singular_values_01110 = np.loadtxt(os.path.join(singular_values_path, f'singular_values_VNN_N3LO_EM500_SLLJT_01110_lambda_1.80_Np_100_np_nocut.dat'))
singular_values_11101 = np.loadtxt(os.path.join(singular_values_path, f'singular_values_VNN_N3LO_EM500_SLLJT_11101_lambda_1.80_Np_100_np_nocut.dat'))
singular_values_11111 = np.loadtxt(os.path.join(singular_values_path, f'singular_values_VNN_N3LO_EM500_SLLJT_11111_lambda_1.80_Np_100_np_nocut.dat'))
singular_values_11121 = np.loadtxt(os.path.join(singular_values_path, f'singular_values_VNN_N3LO_EM500_SLLJT_11121_lambda_1.80_Np_100_np_nocut.dat'))



number_samples = 100

for i in range(number_samples):
    filename_output_np_00001 = f'VNN_N3LO_svPlies{i}_SLLJT_00001_lambda_1.80_Np_100_np_nocut.dat'
    filename_output_np_10010 = f'VNN_N3LO_svPlies{i}_SLLJT_10010_lambda_1.80_Np_100_np_nocut.dat'
    filename_output_np_01110 = f'VNN_N3LO_svPlies{i}_SLLJT_01110_lambda_1.80_Np_100_np_nocut.dat'
    filename_output_np_11101 = f'VNN_N3LO_svPlies{i}_SLLJT_11101_lambda_1.80_Np_100_np_nocut.dat'
    filename_output_np_11111 = f'VNN_N3LO_svPlies{i}_SLLJT_11111_lambda_1.80_Np_100_np_nocut.dat'
    filename_output_np_11121 = f'VNN_N3LO_svPlies{i}_SLLJT_11121_lambda_1.80_Np_100_np_nocut.dat'

    file_output_np_00001 = os.path.join(output_path, filename_output_np_00001)
    file_output_np_10010 = os.path.join(output_path, filename_output_np_10010)
    file_output_np_01110 = os.path.join(output_path, filename_output_np_01110)
    file_output_np_11101 = os.path.join(output_path, filename_output_np_11101)
    file_output_np_11111 = os.path.join(output_path, filename_output_np_11111)
    file_output_np_11121 = os.path.join(output_path, filename_output_np_11121)



    filename_output_pp_00001 = f'VNN_N3LO_svPlies{i}_SLLJT_00001_lambda_1.80_Np_100_pp_nocut.dat'
    #filename_output_pp_10010 = f'VNN_N3LO_svPlies{i}_SLLJT_10010_lambda_1.80_Np_100_pp_nocut.dat'
    #filename_output_pp_01110 = f'VNN_N3LO_svPlies{i}_SLLJT_01110_lambda_1.80_Np_100_pp_nocut.dat'
    filename_output_pp_11101 = f'VNN_N3LO_svPlies{i}_SLLJT_11101_lambda_1.80_Np_100_pp_nocut.dat'
    filename_output_pp_11111 = f'VNN_N3LO_svPlies{i}_SLLJT_11111_lambda_1.80_Np_100_pp_nocut.dat'
    filename_output_pp_11121 = f'VNN_N3LO_svPlies{i}_SLLJT_11121_lambda_1.80_Np_100_pp_nocut.dat'

    file_output_pp_00001 = os.path.join(output_path, filename_output_pp_00001)
    #file_output_pp_10010 = os.path.join(output_path, filename_output_pp_10010)
    #file_output_pp_01110 = os.path.join(output_path, filename_output_pp_01110)
    file_output_pp_11101 = os.path.join(output_path, filename_output_pp_11101)
    file_output_pp_11111 = os.path.join(output_path, filename_output_pp_11111)
    file_output_pp_11121 = os.path.join(output_path, filename_output_pp_11121)



    filename_output_nn_00001 = f'VNN_N3LO_svPlies{i}_SLLJT_00001_lambda_1.80_Np_100_nn_nocut.dat'
    #filename_output_nn_10010 = f'VNN_N3LO_svPlies{i}_SLLJT_10010_lambda_1.80_Np_100_nn_nocut.dat'
    #filename_output_nn_01110 = f'VNN_N3LO_svPlies{i}_SLLJT_01110_lambda_1.80_Np_100_nn_nocut.dat'
    filename_output_nn_11101 = f'VNN_N3LO_svPlies{i}_SLLJT_11101_lambda_1.80_Np_100_nn_nocut.dat'
    filename_output_nn_11111 = f'VNN_N3LO_svPlies{i}_SLLJT_11111_lambda_1.80_Np_100_nn_nocut.dat'
    filename_output_nn_11121 = f'VNN_N3LO_svPlies{i}_SLLJT_11121_lambda_1.80_Np_100_nn_nocut.dat'

    file_output_nn_00001 = os.path.join(output_path, filename_output_nn_00001)
    #file_output_nn_10010 = os.path.join(output_path, filename_output_nn_10010)
    #file_output_nn_01110 = os.path.join(output_path, filename_output_nn_01110)
    file_output_nn_11101 = os.path.join(output_path, filename_output_nn_11101)
    file_output_nn_11111 = os.path.join(output_path, filename_output_nn_11111)
    file_output_nn_11121 = os.path.join(output_path, filename_output_nn_11121)

    x = resampled_LECs[i, :] #current row
    svs_00001 = np.array([x[0], singular_values_00001[1], x[1], singular_values_00001[3], x[2]])
    svs_10010 = np.array([x[3], x[4], x[5], singular_values_10010[3], singular_values_10010[4]])
    svs_01110 = np.array([x[6], x[7], x[8], singular_values_01110[3], singular_values_01110[4]])
    svs_11101 = np.array([x[9], singular_values_11101[1], singular_values_11101[2], x[10], singular_values_11101[4]])  ####CORRECT THIS AND EVERYHING AFTERWARDS
    svs_11111 = np.array([x[11], singular_values_11111[1], x[12], x[13], singular_values_11111[4]])
    svs_11121 = np.array([x[14], x[15], x[16], singular_values_11121[3], singular_values_11121[4]])

    with open(file_output_np_00001, 'w') as f:
        for k in range(5):
            f.write(str(svs_00001[k]))
            f.write('\n')
    with open(file_output_pp_00001, 'w') as f:
        for k in range(5):
            f.write(str(svs_00001[k]))
            f.write('\n')
    with open(file_output_nn_00001, 'w') as f:
        for k in range(5):
            f.write(str(svs_00001[k]))
            f.write('\n')


    with open(file_output_np_10010, 'w') as f:
        for k in range(5):
            f.write(str(svs_10010[k]))
            f.write('\n')

    #with open(file_output_pp_10010, 'w') as f:
    #    for k in range(5):
    #        f.write(str(svs_10010[k]))
    #        f.write('\n')
    #with open(file_output_nn_10010, 'w') as f:
    #    for k in range(5):
    #        f.write(str(svs_10010[k]))
    #        f.write('\n')


    with open(file_output_np_01110, 'w') as f:
        for k in range(5):
            f.write(str(svs_01110[k]))
            f.write('\n')
    #with open(file_output_pp_01110, 'w') as f:
    #    for k in range(5):
    #        f.write(str(svs_01110[k]))
    #        f.write('\n')
    #with open(file_output_nn_01110, 'w') as f:
    #    for k in range(5):
    #        f.write(str(svs_01110[k]))
    #        f.write('\n')


    with open(file_output_np_11101, 'w') as f:
        for k in range(5):
            f.write(str(svs_11101[k]))
            f.write('\n')
    with open(file_output_pp_11101, 'w') as f:
        for k in range(5):
            f.write(str(svs_11101[k]))
            f.write('\n')
    with open(file_output_nn_11101, 'w') as f:
        for k in range(5):
            f.write(str(svs_11101[k]))
            f.write('\n')


    with open(file_output_np_11111, 'w') as f:
        for k in range(5):
            f.write(str(svs_11111[k]))
            f.write('\n')
    with open(file_output_pp_11111, 'w') as f:
        for k in range(5):
            f.write(str(svs_11111[k]))
            f.write('\n')
    with open(file_output_nn_11111, 'w') as f:
        for k in range(5):
            f.write(str(svs_11111[k]))
            f.write('\n')


    with open(file_output_np_11121, 'w') as f:
        for k in range(5):
            f.write(str(svs_11121[k]))
            f.write('\n')
    with open(file_output_pp_11121, 'w') as f:
        for k in range(5):
            f.write(str(svs_11121[k]))
            f.write('\n')
    with open(file_output_nn_11121, 'w') as f:
        for k in range(5):
            f.write(str(svs_11121[k]))
            f.write('\n')