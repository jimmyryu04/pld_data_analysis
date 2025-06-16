import pandas as pd
import os

# Configuration
BASE_DIR = "/qbio/ryu/project/cofolding_data/1E_PRE_analysis/250523"
OUTPUT_DIR = os.path.join(BASE_DIR, "output")
PROTEINS = ["NLucP", "GFP025", "mCherry", "hRLucP", "FLuc2P"]
OUTPUT_FILENAME = "1E_PRE_analysis_summary.txt"

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

output_filepath = os.path.join(OUTPUT_DIR, OUTPUT_FILENAME)

def get_condition_type(name_str):
    """Extracts condition (control, cofolding, masking) from the name string."""
    if name_str.endswith('_CTRL'):
        return 'control'
    elif name_str.endswith('_CFLD'):
        return 'cofolding'
    elif name_str.endswith('_MASK'):
        return 'masking'
    return None

def get_sequence_id(name_str):
    """Extracts the base sequence ID from the name string."""
    if name_str.endswith('_CTRL'):
        return name_str[:-5]  # Remove _CTRL
    elif name_str.endswith('_CFLD'):
        return name_str[:-5]  # Remove _CFLD
    elif name_str.endswith('_MASK'):
        return name_str[:-5]  # Remove _MASK
    return name_str # Should ideally not happen if naming is consistent

all_results = []

for protein in PROTEINS:
    protein_results = [f"Protein: {protein}", "--------------------"]
    filepath = os.path.join(BASE_DIR, f"{protein}_bpp_updated.tsv")

    if not os.path.exists(filepath):
        protein_results.append(f"File not found: {filepath}\n")
        all_results.extend(protein_results)
        all_results.append("====================\n")
        continue

    try:
        df = pd.read_csv(filepath, sep='\t')
    except Exception as e:
        protein_results.append(f"Error reading file {filepath}: {e}\n")
        all_results.extend(protein_results)
        all_results.append("====================\n")
        continue

    # Ensure required columns exist
    required_cols = ['name', 'bp-utp-ms', 'bpp-utp-ms']
    if not all(col in df.columns for col in required_cols):
        protein_results.append(f"Missing one or more required columns ({', '.join(required_cols)}) in {filepath}\n")
        all_results.extend(protein_results)
        all_results.append("====================\n")
        continue
        
    # Convert relevant columns to numeric, coercing errors to NaN
    df['bp-utp-ms'] = pd.to_numeric(df['bp-utp-ms'], errors='coerce')
    df['bpp-utp-ms'] = pd.to_numeric(df['bpp-utp-ms'], errors='coerce')

    # Drop rows where conversion to numeric failed for key metrics
    df.dropna(subset=['bp-utp-ms', 'bpp-utp-ms'], inplace=True)

    df['condition'] = df['name'].apply(get_condition_type)
    df['sequence_id'] = df['name'].apply(get_sequence_id)

    # Filter out rows where condition or sequence_id could not be determined properly
    df = df[df['condition'].notna()]
    df = df[~df['sequence_id'].str.endswith(('_CTRL', '_CFLD', '_MASK'))] # Ensure sequence_id is clean

    protein_results.append("Averages per condition:")
    conditions_data = {}
    for condition_name in ['control', 'cofolding', 'masking']:
        subset = df[df['condition'] == condition_name]
        conditions_data[condition_name] = subset[['sequence_id', 'bp-utp-ms', 'bpp-utp-ms']].set_index('sequence_id')
        
        if not subset.empty:
            count = len(subset)
            avg_bp_utp_ms = subset['bp-utp-ms'].mean()
            avg_bpp_utp_ms = subset['bpp-utp-ms'].mean()
            protein_results.append(f"  Condition: {condition_name}")
            protein_results.append(f"    Number of sequences: {count}")
            protein_results.append(f"    Average bp-utp-ms: {avg_bp_utp_ms:.4f}")
            protein_results.append(f"    Average bpp-utp-ms: {avg_bpp_utp_ms:.4f}")
        else:
            protein_results.append(f"  Condition: {condition_name}")
            protein_results.append(f"    Number of sequences: 0")
            protein_results.append(f"    Average bp-utp-ms: N/A")
            protein_results.append(f"    Average bpp-utp-ms: N/A")
    protein_results.append("") # Newline for spacing

    # Comparisons
    control_df = conditions_data.get('control')
    masking_df = conditions_data.get('masking')
    cofolding_df = conditions_data.get('cofolding')

    # Masking vs. Control
    protein_results.append("Masking vs. Control Comparison:")
    if control_df is not None and not control_df.empty and masking_df is not None and not masking_df.empty:
        merged_ctrl_mask = control_df.join(masking_df, lsuffix='_ctrl', rsuffix='_mask', how='inner')
        if not merged_ctrl_mask.empty:
            total_compared = len(merged_ctrl_mask)
            # Ensure columns from suffixes exist after join
            if 'bp-utp-ms_ctrl' in merged_ctrl_mask.columns and 'bp-utp-ms_mask' in merged_ctrl_mask.columns:
                mask_gt_ctrl_bp = (merged_ctrl_mask['bp-utp-ms_mask'] > merged_ctrl_mask['bp-utp-ms_ctrl']).sum()
                equal_bp_count = (merged_ctrl_mask['bp-utp-ms_mask'] == merged_ctrl_mask['bp-utp-ms_ctrl']).sum()
                ratio_bp_str = f"{mask_gt_ctrl_bp}/{total_compared} ({mask_gt_ctrl_bp/total_compared:.4f})" if total_compared > 0 else "0/0 (N/A)"
                ratio_equal_bp_str = f"{equal_bp_count}/{total_compared} ({equal_bp_count/total_compared:.4f})"
            else: # Handle case where one of the columns might be all NaN and dropped or not present
                ratio_bp_str = "N/A (Required columns missing for comparison)"
                ratio_equal_bp_str = "N/A (Required columns missing for comparison)"

            if 'bpp-utp-ms_ctrl' in merged_ctrl_mask.columns and 'bpp-utp-ms_mask' in merged_ctrl_mask.columns:
                mask_gt_ctrl_bpp = (merged_ctrl_mask['bpp-utp-ms_mask'] > merged_ctrl_mask['bpp-utp-ms_ctrl']).sum()
                equal_bpp_count = (merged_ctrl_mask['bpp-utp-ms_mask'] == merged_ctrl_mask['bpp-utp-ms_ctrl']).sum()
                ratio_bpp_str = f"{mask_gt_ctrl_bpp}/{total_compared} ({mask_gt_ctrl_bpp/total_compared:.4f})" if total_compared > 0 else "0/0 (N/A)"
                ratio_equal_bpp_str = f"{equal_bpp_count}/{total_compared} ({equal_bpp_count/total_compared:.4f})"
            else:
                ratio_bpp_str = "N/A (Required columns missing for comparison)"
                ratio_equal_bpp_str = "N/A (Required columns missing for comparison)"


            protein_results.append(f"  Total sequences compared: {total_compared}")
            protein_results.append(f"  Masking > Control (bp-utp-ms): {ratio_bp_str}")
            protein_results.append(f"  Masking = Control (bp-utp-ms): {ratio_equal_bp_str}")
            protein_results.append(f"  Masking > Control (bpp-utp-ms): {ratio_bpp_str}")
            protein_results.append(f"  Masking = Control (bpp-utp-ms): {ratio_equal_bpp_str}")
        else:
            protein_results.append("  No common sequences found for comparison.")
    else:
        protein_results.append("  Insufficient data for Control or Masking condition.")
    protein_results.append("") # Newline for spacing

    
    all_results.extend(protein_results)
    all_results.append("\n====================\n")


# Write results to file
with open(output_filepath, 'w') as f_out:
    for line in all_results:
        f_out.write(line + '\n')

print(f"Analysis complete. Results saved to: {output_filepath}")