#!/bin/bash
#SBATCH --job-name=A71E50_bpp
#SBATCH --output=/home/ryu/project/ParalinearDesign/log/slurm-%x-%j.out
#SBATCH --error=/home/ryu/project/ParalinearDesign/log/slurm-%x-%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=55
#SBATCH --time=2-00:00:00
#SBATCH --partition=cpu

source /opt/ohpc/pub/anaconda3/etc/profile.d/conda.sh
conda activate pld

DIR="/home/ryu/project/ParalinearDesign/250611_50_cds-A7_1E"
COF="/home/ryu/project/ParalinearDesign/codon-opt-factors"
PROTEINS="NLucP"

# NLucP GFP025 mCherry hRLucP FLuc2P
                    # --with-mask \
                    #  --cds-only \
                    #  --process 1 \
                    #  --threads 50 \
                    #  --num-random-samples 50 \
                    #  --random-rejection-rate 0.02

for protein in ${PROTEINS}; do
  # python3 /home/ryu/project/ParalinearDesign/pld-tools/src/run_pld_random_sampling.py --utr3-fwd ${DIR}/pre.fasta \
  #                    --utr3-bwd ${DIR}/utr_bwd.fasta \
  #                    --cds-fasta ${DIR}/proteins.faa \
  #                    --output-dir ${DIR} \
  #                    --lambda 2,3,4 \
  #                    --cds-only  \
  #                    --process 1 \
  #                    --threads 20 \
  #                    --num-random-samples 500 \
  #                    --random-rejection-rate 0.02

  python3 /home/ryu/project/ParalinearDesign/pld-tools/src/parse_pld_result_random_sampling.py --input-utp "${DIR}/${protein}.fasta" \
                              --output-dir "${DIR}"

  python3 /home/ryu/project/ParalinearDesign/pld-tools/src/analyze_bp.py --input-tsv "${DIR}/${protein}.tsv"\
                        --output-dir "${DIR}"

  python3 /home/ryu/project/ParalinearDesign/pld-tools/src/analyze_bpp.py --input-tsv "${DIR}/${protein}.tsv"\
                         --bpp-dir ${DIR}/bpp \
                         --output-dir "${DIR}" \
                         --len-5utr 53

  # python3 /home/ryu/project/ParalinearDesign/pld-tools/src/merge_cof_metrics.py --input-utp-tsv ${DIR}/${protein}_bpp.tsv\
  #                              --output-dir ${DIR} \
  #                              --output-prefix "${protein}" \
  #                              --cof-repo ${COF} \
  #                              --merge-metrics aup_unweighted\
  #                              --merge-metrics aup_U3A1o5\
  #                              --merge-metrics degscore\
  #                              --merge-metrics start_str\
  #                              --merge-metrics start_app

  # python3 /home/ryu/project/ParalinearDesign/pld-tools/src/filter.py --input-tsv "${COF}/${protein}_merged.tsv"\
  #                   --output-dir "${DIR}/"
done
