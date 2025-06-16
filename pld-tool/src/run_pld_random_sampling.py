# pld-tools/src/run_pld_random_sampling.py

import click
import concurrent.futures
import numpy as np
import pandas as pd
import subprocess
import RNA
from Bio import SeqIO
from collections import defaultdict
from itertools import product
from pathlib import Path
from rich.console import Console
from rich.progress import (
    Progress, BarColumn, TextColumn,
    MofNCompleteColumn, TimeElapsedColumn,
    TimeRemainingColumn
)

console = Console()

HBB5 = "AGGACAUUUGCUUCUGACACAACUGUGUUCACUAGCAACCUCAAACAGACACC"
CODON_USAGE_FILE = "/home/ryu/project/ParalinearDesign/pld-tools/resources/codon_usage_freq_table_human.csv"

def load_query_data(fasta_path: Path) -> tuple[dict, dict]:
    """Load query sequences and structures from extended FASTA format."""
    query_seq = {}
    query_struct = {}
    with open(fasta_path) as f:
        lines = [line.strip() for line in f if line.strip()]
        for i in range(0, len(lines), 3):
            name = lines[i][1:]
            seq = lines[i + 1].upper().replace("T", "U")
            struct = lines[i + 2]
            query_seq[name] = seq
            query_struct[name] = struct
    return query_seq, query_struct

def calculate_cai(seq, codonusage_filename):
    """Calculate Codon Adaptation Index (CAI) from a given CDS using a codon usage table."""
    freqtbl = pd.read_csv(codonusage_filename, names=['codon', 'aa', 'freq'], comment='#')
    freqtbl['freq'] = freqtbl['freq'].replace(0, 1e-100)
    freqtbl['cai'] = freqtbl['freq'] / freqtbl.groupby('aa')['freq'].transform('max')
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    weights = [freqtbl.loc[freqtbl['codon'] == codon, 'cai'].values[0]
               for codon in codons if codon in freqtbl['codon'].values]
    return round(np.prod(weights) ** (1 / len(weights)), 3) if weights else 0.0

def parse_pldesign_fasta_output(fasta_output):
    """Parse FASTA-formatted pldesign output, potentially with multiple entries."""
    sequences_data = []
    lines = fasta_output.strip().splitlines()
    if not lines:
        return sequences_data
        
    for i in range(0, len(lines), 2):
        if i + 1 < len(lines):
            header = lines[i]
            sequence = lines[i+1]
            if not header.startswith('>'): continue
            sequences_data.append(sequence.strip())
            
    return sequences_data

def run_pldesign(params):
    """Run the pldesign command with given parameters and return extracted sequence."""
    result = subprocess.run(["pldesign"] + params, capture_output=True, text=True)
    return parse_pldesign_fasta_output(result.stdout)

def process_case(case_name, cds_seq, utr_seq):
    """Construct full mRNA, fold it, and compute MFE and CAI."""
    full_seq = HBB5.lower() + cds_seq + (utr_seq.lower() if utr_seq else "")
    struct, mfe = RNA.fold(full_seq)
    cai = calculate_cai(cds_seq, CODON_USAGE_FILE)
    return case_name, mfe, cai, full_seq, struct

def get_mask_region(prot_seq, utr_id, utr_seq):
    """Find motif match region in UTR and return masking region strings, expanded by 3 nts."""
    query_seq, query_struct = load_query_data(Path("/home/ryu/project/ParalinearDesign/pld-tools/resources/motif_data.fasta"))

    for name, qseq in query_seq.items():
        if name in utr_id:
            q_upper = qseq.upper().replace("T", "U")
            utr_upper = utr_seq.upper().replace("T", "U")
            
            motif_start_in_utr_0idx = utr_upper.find(q_upper)
            
            if motif_start_in_utr_0idx != -1:
                cds_len = len(prot_seq) * 3
                original_motif_start_0idx = motif_start_in_utr_0idx + cds_len
                original_motif_end_0idx = original_motif_start_0idx + len(qseq) - 1
                
                sequence_full_len = cds_len + len(utr_seq)

                expanded_mask_start_0idx = original_motif_start_0idx - 3
                expanded_mask_end_0idx = original_motif_end_0idx + 3

                clamped_expanded_mask_start_0idx = max(0, expanded_mask_start_0idx)
                clamped_expanded_mask_end_0idx = min(sequence_full_len - 1, expanded_mask_end_0idx)

                return (
                    f"{clamped_expanded_mask_start_0idx},{clamped_expanded_mask_end_0idx},1,{original_motif_start_0idx}",
                    f"{clamped_expanded_mask_start_0idx},{clamped_expanded_mask_end_0idx},{original_motif_end_0idx},{sequence_full_len}",
                    f"1,{original_motif_start_0idx},{clamped_expanded_mask_start_0idx},{clamped_expanded_mask_end_0idx}",
                    f"{original_motif_end_0idx},{sequence_full_len},{clamped_expanded_mask_start_0idx},{clamped_expanded_mask_end_0idx}"
                )
    return None

def run_case(args):
    """Run one case of pldesign, process result, and append it to result list and tmp file."""
    case_type, utr_id1, utr_id2, utr_seq, lam, prot_id, prot_seq, result_list, progress, task_id, tmp_path, pld_sampling_params = args
    
    base_cmd = ["-p", prot_seq, "-l", str(lam), "--fasta", "--exclude-start", "5"]
    
    if pld_sampling_params['num_random_samples'] > 1:
        base_cmd.extend([
            "--num-random-samples", str(pld_sampling_params['num_random_samples']),
            "--random-rejection-rate", str(pld_sampling_params['random_rejection_rate']),
            "--threads", str(pld_sampling_params['threads'])
        ])

    pldesign_results = []
    if case_type == 'cds':
        case_name = f"{prot_id}_{utr_id1}_{utr_id2}_lam{lam}_CTRL"
        pldesign_results = run_pldesign(base_cmd)
    elif case_type == 'utr':
        case_name = f"{prot_id}_{utr_id1}_{utr_id2}_lam{lam}_CFLD"
        cmd = base_cmd + ["-u", utr_seq]
        pldesign_results = run_pldesign(cmd)
    elif case_type == 'mask':
        case_name = f"{prot_id}_{utr_id1}_{utr_id2}_lam{lam}_MASK"
        mask_regions = get_mask_region(prot_seq, utr_id1, utr_seq)
        if mask_regions:
            cmd = base_cmd + ["-u", utr_seq]
            for region in mask_regions:
                cmd += ["-m", region]
            pldesign_results = run_pldesign(cmd)
        else:
            progress.update(task_id, advance=1, refresh=True)
            return
    else:
        progress.update(task_id, advance=1, refresh=True)
        return
        
    for i, seq_res in enumerate(pldesign_results):
        if case_type in ['utr', 'mask']:
            cds_seq = seq_res[:-len(utr_seq)]
        else: 
            cds_seq = seq_res

        case_name_with_sample = f"{case_name}_sample{i+1}"
        
        c_name, mfe, cai, full_seq, struct = process_case(case_name_with_sample, cds_seq, utr_seq)

        result_list.append((prot_id, case_type, c_name, mfe, cai, full_seq, struct))

        with open(tmp_path, "a") as tf:
            tf.write(f"{prot_id}\t{case_type}\t{c_name}\t{len(full_seq)}\t{mfe:.2f}\t{cai:.3f}\t{full_seq}\t{struct}\n")

    progress.update(task_id, advance=1, refresh=True)

def write_sorted_outputs(output_prefix, result_list):
    """Sort and write results to FASTA-formatted files grouped by protein ID."""
    result_list.sort(key=lambda x: (x[0], x[2]))

    results_by_prot = defaultdict(list)
    for prot_id, case_type, name, mfe, cai, seq, struct in result_list:
        results_by_prot[prot_id].append((name, mfe, cai, seq, struct))

    for prot_id, entries in results_by_prot.items():
        fasta_path = f"{output_prefix}{prot_id}.fasta"
        with open(fasta_path, "w") as ff:
            for name, mfe, cai, seq, struct in entries:
                ff.write(f">{name}\t{mfe:.2f}\t{cai:.3f}\n{seq}\n{struct}\n")
        console.print(f"[bold white]=> FASTA saved: {fasta_path}")

def prepare_task_args(
    utr1_records,
    utr2_records,
    lambda_list,
    selected_cases,
    prot_id,
    prot_seq,
    progress,
    tmp_path,
    pld_sampling_params
) -> list[tuple]:
    """Prepare argument tuples for parallel pldesign execution for a single protein."""
    task_args = []
    task_total = len(utr1_records) * len(utr2_records) * len(lambda_list) * len(selected_cases)
    task_id = progress.add_task(f"{prot_id:>10}...", total=task_total)

    for record1, record2 in product(utr1_records, utr2_records):
        utr_seq = str(record1.seq) + str(record2.seq)
        for lam in lambda_list:
            for case_type in selected_cases:
                task_args.append((
                    case_type, record1.id, record2.id, utr_seq, lam,
                    prot_id, prot_seq, None,
                    progress, task_id, tmp_path, pld_sampling_params
                ))
    return task_args

@click.command()
@click.option('--utr3-fwd', type=click.Path(exists=True), required=True,
              help='FASTA files with 3\' UTR candidates.')
@click.option('--utr3-bwd', type=click.Path(exists=True), required=True,
              help='FASTA files with 3\' UTR candidates.')
@click.option('--cds-fasta', type=click.Path(exists=True), required=True,
              help='FASTA file with protein sequences for CDS generation.')
@click.option('--output-dir', type=str, default=".",
              help='Path of output directory.')
@click.option('--output-prefix', type=str, default=None,
              help='Prefix for output filenames.')
@click.option('--cds-only', is_flag=True, help='Run CDS-only optimization.')
@click.option('--with-utr', is_flag=True, help='Run optimization with UTR.')
@click.option('--with-mask', is_flag=True, help='Run optimization with UTR and masking.')
@click.option('--lambda', 'lambda_values', type=str, default="2,3,4",
              help='Comma-separated lambda values (e.g., "2,3,4").')
@click.option('--num-random-samples', '-r', type=int, default=1,
              help='Number of random samples to generate.')
@click.option('--random-rejection-rate', type=float, default=0.01,
              help='Rejection rate for random sampling.')
@click.option('--threads', '-t', type=int, default=4, help='Number of threads to use for pldesign.')
@click.option('--process', '-p', type=int, default=4, help='Number of parallel processes to run this script.')
def main(utr3_fwd, utr3_bwd, cds_fasta, output_dir, output_prefix,
         cds_only, with_utr, with_mask, lambda_values, num_random_samples,
         random_rejection_rate, threads, process):
    lambda_list = [int(l.strip()) for l in lambda_values.split(",") if l.strip().isdigit()]

    pld_sampling_params = {
        'num_random_samples': num_random_samples,
        'random_rejection_rate': random_rejection_rate,
        'threads': threads,
    }

    utr1_records = list(SeqIO.parse(utr3_fwd, "fasta"))
    utr2_records = list(SeqIO.parse(utr3_bwd, "fasta"))
    cds_records = list(SeqIO.parse(cds_fasta, "fasta"))

    selected_cases = set()
    if not any([cds_only, with_utr, with_mask]):
        selected_cases = {'cds', 'utr', 'mask'}
    else:
        if cds_only:
            selected_cases.add('cds')
        if with_utr:
            selected_cases.add('utr')
        if with_mask:
            selected_cases.add('mask')

    if output_prefix:
        output_prefix = output_prefix + "_"
    else:
        output_prefix = ""
    tmp_path = f"{output_dir}/tmp.txt"
    with open(tmp_path, "w") as tf:
        tf.write("prot_id\tcase_type\tcase_name\tlength\tmfe\tcai\tsequence\tstructure\n")

    results = []

    with Progress(
        TextColumn("{task.description}"),
        BarColumn(bar_width=None, complete_style="bold green"),
        MofNCompleteColumn(),
        TextColumn("•"),
        TimeElapsedColumn(),
        TextColumn("•"),
        TimeRemainingColumn(),
        console=console,
    ) as progress:
        for cds_record in cds_records:
            prot_id = cds_record.id
            prot_seq = str(cds_record.seq)

            task_args = prepare_task_args(
                utr1_records, utr2_records, lambda_list,
                selected_cases, prot_id, prot_seq,
                progress, tmp_path, pld_sampling_params
            )
            for i in range(len(task_args)):
                task_args[i] = list(task_args[i])
                task_args[i][7] = results

            with concurrent.futures.ThreadPoolExecutor(max_workers=process) as executor:
                futures = [executor.submit(run_case, tuple(arg)) for arg in task_args]
                for future in concurrent.futures.as_completed(futures):
                    try:
                        future.result()
                    except Exception as e:
                        console.print(f"[red]Exception in run_case: {e}[/red]")

    write_sorted_outputs(f"{output_dir}/{output_prefix}", results)

if __name__ == '__main__':
    main()