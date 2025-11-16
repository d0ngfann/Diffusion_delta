# SLURM 작업 배열로 시뮬레이션을 병렬 실행하기 위한 배치 스크립트
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --job-name=cf_fill1
#SBATCH --mem=100G
#SBATCH --partition=short
#SBATCH --output=//home/wan.a/complex_contagion_repo/output_job/slurm-%A_%a.out

module load anaconda3/2022.01

source activate centola_riff_env
#echo $python 3.7.7

declare -a commands
# replace batch_main_perc0.conf with any other configuration file
readarray -t commands < batch_main_perc0.conf # Exclude newline.
eval ${commands[$SLURM_ARRAY_TASK_ID]}


# to run in console: sbatch --array=0-[ARRAY LENGTH] run_seeding.sh, 0-611
# to see job status: squeue -u wan.a
# to cancel job: scancel job_number