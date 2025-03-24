import numpy as np
import os
import random
import matplotlib.pyplot as plt
file_path = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/samples'


number_walkers = 100
number_steps = 10000
burn_in = 5000

number_parameters = 19

output_path = '/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/subsamples'


#samples_filename = 'MCMC_EKM_EM_3N_10000_no_cis_25MeV_coupled_A.txt'
samples_filename = 'MCMC_EKM_EM_3N_10000_no_cis_coupled_1MeV.txt'
samples_file = os.path.join(file_path, samples_filename)
samples = np.loadtxt(samples_file)

#phaseshifts_filename = 'EKM_EM_phaseshift_1_optimized_17_params_50MeV.dat'
#phaseshifts_filename = 'EKM_EM_phaseshift_3N_new_likelihood_fixed.dat'
#phaseshifts_file = os.path.join(file_path, phaseshifts_filename)
#phaseshifts = np.loadtxt(phaseshifts_file)[number_walkers * 6:, :]





print('shape samples: ', np.shape(samples))

separate_arrays_samples = []

for i in range(number_walkers):

    sliced_array_samples = samples[i::number_walkers][burn_in:]
    separate_arrays_samples.append(sliced_array_samples)

separate_arrays_samples = np.array(separate_arrays_samples)


print(np.shape(separate_arrays_samples))

resampled_samples = []

for k in range(number_walkers):

    m = random.randint(0, number_steps - burn_in -1) # random number m in range number_steps - burn_in
    print(m)
    print(k)
    resampled_samples.append(separate_arrays_samples[k, m, :])



#write resampled stuff out to file


resampled_samples = np.array(resampled_samples)
print(np.shape(resampled_samples))


output_filename = '3N_10000_no_cis_1MeV_coupled.dat'
output_file = os.path.join(output_path, output_filename)
with open(output_file, 'w') as f:
    for individual_sample in range(number_walkers):
        for parameter in range(number_parameters):
            f.write(str(resampled_samples[individual_sample, parameter]))
            f.write(' ')
        f.write('\n')



plt.scatter(resampled_samples[:, 17], resampled_samples[:, 18])
#plt.savefig('resampled_parameters_8000_no_cis_NEW.pdf')
plt.show()


dont run !!

