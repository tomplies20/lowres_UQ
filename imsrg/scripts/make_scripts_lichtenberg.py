import numpy as np

likelihood = '25MeV'

def create_bash_script(
        job_name,
        NN_filename,
        path,
        hw,
        A,
        reference,
        valence_space,
        BetaCM,
        hwBetaCM,
        file2e1max,
        file2e2max,
        file2lmax,
        file3e1max,
        file3e2max,
        file3e3max,
        THREE_bme_type,
        no2b_precision,
        c1,
        c3,
        c4,
        cD,
        cE,
        basis,
        use_HF_reference_in_NAT,
        emax,
        e3max,
        emax_imsrg,
        method,
        smax,
        dsmax,
        ds_0,
        omega_norm_max,
        eta_criterion,
        write_omega,
        Operators,
        number_job,
):

    job_name = f'./script_files_lichtenberg_{reference}/hw_{hw}_emax_{emax}_A_{A}_reference_{reference}_likelihood_{likelihood}_number_job_{number_job}_emax_imsrg_10'
    #_c1_{c1}_c3_{c3}_c4_{c4}_cD_{cD}_cE_{cE}
    # Define the content of the bash script
    bash_script_content = f"""#!/bin/bash
#SBATCH -A project02474
#SBATCH --job-name={job_name[27:]}
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3800
#SBATCH --time=8:00:00
##SBATCH --mail-user=tplies@theorie.ikp.physik.tu-darmstadt.de
##SBATCH --mail-type=END,FAIL
#SBATCH --output={path}/output/{job_name}_STDOUT.txt
#SBATCH --error={path}/output/{job_name}_STDERR.txt

module --ignore-cache load intel boost gsl cmake
module refresh
export OMP_NUM_THREADS=48 
cd {path} 
mkdir -p /work/scratch/tp41hatu/data/scratch/scratch_{job_name} 
mkdir -p {path}/output/output_{job_name} 
{path}/src/build/imsrg++_3n_combine \\
scratch=/work/scratch/tp41hatu/data/scratch/scratch_{job_name}  \\
hw={hw} A={A} reference={reference} valence_space={valence_space}   \\
BetaCM={BetaCM} hwBetaCM={hwBetaCM}  \\
2bme=/work/scratch/tp41hatu/NN_matrix_elements/{NN_filename}  \\
file2e1max={file2e1max} file2e2max={file2e2max} file2lmax={file2lmax}  \\
3bme_c1=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_1.0_c3_0.0_c4_0.0_cD_0.0_cE_0.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \\
3bme_c3=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_0.0_c3_1.0_c4_0.0_cD_0.0_cE_0.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \\
3bme_c4=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_0.0_c3_0.0_c4_1.0_cD_0.0_cE_0.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \\
3bme_cD=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_0.0_c3_0.0_c4_0.0_cD_1.0_cE_0.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \\
3bme_cE=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_0.0_c3_0.0_c4_0.0_cD_0.0_cE_1.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \\
file3e1max={file3e1max} file3e2max={file3e2max} file3e3max={file3e3max} 3bme_type={THREE_bme_type} no2b_precision={no2b_precision}  \\
c1={c1} c3={c3} c4={c4} cD={cD} cE={cE}  \\
basis={basis} use_HF_reference_in_NAT={use_HF_reference_in_NAT} emax={emax}  e3max={e3max} \\
emax_imsrg={emax_imsrg} \\
method={method} smax={smax} dsmax={dsmax} ds_0={ds_0} omega_norm_max={omega_norm_max} eta_criterion={eta_criterion}  \\
write_omega={write_omega} \\
Operators={Operators} \\
flowfile={path}/output/{job_name}_FLOWFILE.txt \\
intfile={path}/output/output_{job_name}_INTFILE 


rm -r /work/scratch/tp41hatu/data/scratch/scratch_{job_name} 

"""


    with open(job_name + '.sh', "w") as file:
        file.write(bash_script_content)

    print(f"Bash script '{job_name}' created successfully!")



number_samples = 100

samples_path = f'/Users/pleazy/PycharmProjects/uncertainty_quantification/library/results/subsamples/3N_10000_no_cis_{likelihood}_coupled.dat'
samples = np.loadtxt(samples_path)
for i in range(number_samples):

    #c1 = samples[i, 17]
    #c3 = samples[i, 18]
    #c4 = samples[i, 19]
    #cD = samples[i, 20]
    #cE = samples[i, 21]
    c1 = -0.81
    c3 = -3.2
    c4 = 5.4
    cD = samples[i, 17]
    cE = samples[i, 18]
    #cD = 1.264
    #cE = -0.12
    print(c1)
    create_bash_script(
        NN_filename = f'VNN_N3LO_Plies{i}_lambda_1.80_hw16.00_emax_14.me2j ', #TwBME-HO_NN-only_N3LO_EM500_srg1.8_hw16.0_emax18_e2max36.me2j.gz
        path = '/home/tp41hatu/code/imsrg_3N/imsrg',
        hw = 16,
        A = 16,
        reference = 'O16',
        valence_space ='O16',
        BetaCM = 0.0,
        hwBetaCM = 16,
        file2e1max = 18,
        file2e2max = 36,
        file2lmax = 18,
        file3e1max = 18,
        file3e2max = 36,
        file3e3max = 24,
        THREE_bme_type = 'no2b',
        no2b_precision = 'single',
        c1 = c1,
        c3 = c3,
        c4 = c4,
        cD =cD,
        cE = cE,
        basis = 'NAT',
        use_HF_reference_in_NAT = 'true',
        emax = 14,
        e3max = 24,
        emax_imsrg = 10,
        method = 'magnus',
        smax = 100.0,
        dsmax = 0.5,
        ds_0 = 0.01,
        omega_norm_max = 0.25,
        eta_criterion = 0.001,
        write_omega = 'true',
        Operators = 'Rp2c,Rn2c,Rso2',
        job_name='',
        number_job=i,
    )