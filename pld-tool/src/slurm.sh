#!/bin/bash
#SBATCH --job-name=merge_fasta
#SBATCH --output=/home/ryu/project/ParalinearDesign/log/slurm-%x-%j.out
#SBATCH --error=/home/ryu/project/ParalinearDesign/log/slurm-%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=24:00:00
#SBATCH --partition=cpu

# 로그 디렉토리 생성
mkdir -p log

# Conda 환경 활성화
echo "Activating Conda environment..."
source /opt/ohpc/pub/anaconda3/etc/profile.d/conda.sh
conda activate pld

# # 파일 경로 변수 설정
# # 실제 파일이 있는 경로로 수정하세요.
# BASE_DIR="/home/ryu/project/ParalinearDesign/250612_500_merged" # 예시 경로
# CDS_FASTA="${BASE_DIR}/250610_50_cds-k3,4,5/NLucP.fasta"
# MASK_FASTA="${BASE_DIR}/250610_test/NLucP.fasta"
# OUTPUT_FASTA="${BASE_DIR}/250610_test/NLucP_merged.fasta"

# Python 스크립트 실행
echo "Starting FASTA merge and sort..."
python3 /home/ryu/project/ParalinearDesign/pld-tools/src/merge_fasta.py \
    --cds-file "/home/ryu/project/ParalinearDesign/250612_500_cds_lam3_4_k5/NLucP.fasta" \
    --mask-file "/home/ryu/project/ParalinearDesign/250614_500_lam_all_merged_k5/NLucP.fasta" \
    --output-file "/home/ryu/project/ParalinearDesign/250614_500_lam_all_merged_k5/NLucP.fasta"

echo "Job finished successfully."
