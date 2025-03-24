#!/bin/bash
#SBATCH -A project02474
#SBATCH --job-name=O24/hw_16_emax_14_A_24_reference_O24_likelihood_25MeV_number_job_27_emax_imsrg_10
#SBATCH --nodes=1
#SBATCH -n 1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3800
#SBATCH --time=8:00:00
##SBATCH --mail-user=tplies@theorie.ikp.physik.tu-darmstadt.de
##SBATCH --mail-type=END,FAIL
#SBATCH --output=/home/tp41hatu/code/imsrg_3N/imsrg/output/./script_files_lichtenberg_O24/hw_16_emax_14_A_24_reference_O24_likelihood_25MeV_number_job_27_emax_imsrg_10_STDOUT.txt
#SBATCH --error=/home/tp41hatu/code/imsrg_3N/imsrg/output/./script_files_lichtenberg_O24/hw_16_emax_14_A_24_reference_O24_likelihood_25MeV_number_job_27_emax_imsrg_10_STDERR.txt

module load intel boost gsl cmake
export OMP_NUM_THREADS=48 
cd /home/tp41hatu/code/imsrg_3N/imsrg 
mkdir -p /work/scratch/tp41hatu/data/scratch/scratch_./script_files_lichtenberg_O24/hw_16_emax_14_A_24_reference_O24_likelihood_25MeV_number_job_27_emax_imsrg_10 
mkdir -p /home/tp41hatu/code/imsrg_3N/imsrg/output/output_./script_files_lichtenberg_O24/hw_16_emax_14_A_24_reference_O24_likelihood_25MeV_number_job_27_emax_imsrg_10 
/home/tp41hatu/code/imsrg_3N/imsrg/src/build/imsrg++_3n_combine \
scratch=/work/scratch/tp41hatu/data/scratch/scratch_./script_files_lichtenberg_O24/hw_16_emax_14_A_24_reference_O24_likelihood_25MeV_number_job_27_emax_imsrg_10  \
hw=16 A=24 reference=O24 valence_space=O24   \
BetaCM=0.0 hwBetaCM=16  \
2bme=/work/scratch/tp41hatu/NN_matrix_elements/VNN_N3LO_Plies27_lambda_1.80_hw16.00_emax_14.me2j   \
file2e1max=18 file2e2max=36 file2lmax=18  \
3bme_c1=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_1.0_c3_0.0_c4_0.0_cD_0.0_cE_0.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \
3bme_c3=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_0.0_c3_1.0_c4_0.0_cD_0.0_cE_0.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \
3bme_c4=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_0.0_c3_0.0_c4_1.0_cD_0.0_cE_0.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \
3bme_cD=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_0.0_c3_0.0_c4_0.0_cD_1.0_cE_0.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \
3bme_cE=/work/scratch/tp41hatu/data/NO2B_ThBME_3NFJmax15_c1_0.0_c3_0.0_c4_0.0_cD_0.0_cE_1.0_NonLocal4_394.6539576_IS_hw16.0_ms18_36_24.stream.bin   \
file3e1max=18 file3e2max=36 file3e3max=24 3bme_type=no2b no2b_precision=single  \
c1=-0.81 c3=-3.2 c4=5.4 cD=1.5906690909838397 cE=-0.11917058880909742  \
basis=NAT use_HF_reference_in_NAT=true emax=14  e3max=24 \
emax_imsrg=10 \
method=magnus smax=100.0 dsmax=0.5 ds_0=0.01 omega_norm_max=0.25 eta_criterion=0.001  \
write_omega=true \
Operators=Rp2c,Rn2c,Rso2 \
flowfile=/home/tp41hatu/code/imsrg_3N/imsrg/output/./script_files_lichtenberg_O24/hw_16_emax_14_A_24_reference_O24_likelihood_25MeV_number_job_27_emax_imsrg_10_FLOWFILE.txt \
intfile=/home/tp41hatu/code/imsrg_3N/imsrg/output/output_./script_files_lichtenberg_O24/hw_16_emax_14_A_24_reference_O24_likelihood_25MeV_number_job_27_emax_imsrg_10_INTFILE 


rm -r /work/scratch/tp41hatu/data/scratch/scratch_./script_files_lichtenberg_O24/hw_16_emax_14_A_24_reference_O24_likelihood_25MeV_number_job_27_emax_imsrg_10 

