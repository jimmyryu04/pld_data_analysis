import os

def parse_fasta_file(filepath):
    """
    Parses a FASTA file where each entry is exactly 3 lines:
    >Header
    Sequence
    Structure
    Returns a list of tuples, each tuple being (header, sequence, structure).
    Skips empty lines.
    """
    entries = []
    try:
        with open(filepath, 'r') as infile:
            lines = [line.strip() for line in infile if line.strip()] # Read all non-empty stripped lines
        
        if not lines:
            print(f"Warning: File {filepath} is empty.")
            return entries

        i = 0
        while i < len(lines):
            if not lines[i].startswith('>'):
                print(f"Warning: Expected header starting with '>' but got '{lines[i]}' in {filepath} at original line index approx {i}. Skipping this supposed entry.")
                # Try to find next header
                j = i + 1
                while j < len(lines) and not lines[j].startswith('>'):
                    j += 1
                i = j
                continue

            header = lines[i]
            
            if i + 2 < len(lines):
                sequence = lines[i+1]
                structure = lines[i+2]
                # Basic check if sequence or structure accidentally look like a header
                if sequence.startswith('>') or structure.startswith('>'):
                    print(f"Warning: Malformed entry for header '{header}' in {filepath}. Sequence or structure line starts with '>'. Skipping.")
                    i += 1 # Advance by one to allow resync
                    continue
                entries.append((header, sequence, structure))
                i += 3
            else:
                print(f"Warning: Incomplete entry for header '{header}' in {filepath}. Skipping.")
                i += 1 # Advance by one to allow resync for the last entry if incomplete
                break # Or i += 1 if there could be more malformed lines

    except FileNotFoundError:
        print(f"Error: File not found {filepath}")
        return None # Indicate error
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None # Indicate error
    return entries

def main():
    old_fasta_path = "/home/ryu/project/ParalinearDesign/250523/FLuc2P.fasta"
    new_mask_source_path = "/home/ryu/project/ParalinearDesign/250526-change-motif-sites-3/FLuc2P.fasta"
    
    # 최종 결과 파일 경로 설정 (예: ParalinearDesign 폴더 내)
    output_fasta_path = "/home/ryu/project/ParalinearDesign/FLuc2P_updated.fasta"

    print(f"Loading non-MASK sequences from: {old_fasta_path}")
    old_entries = parse_fasta_file(old_fasta_path)
    if old_entries is None:
        return # Stop if file reading failed

    print(f"Loading new MASK sequences from: {new_mask_source_path}")
    new_mask_entries = parse_fasta_file(new_mask_source_path)
    if new_mask_entries is None:
        return # Stop if file reading failed

    final_entries = []

    # 1. old_fasta_path에서 _MASK가 아닌 서열들만 선택
    for header, sequence, structure in old_entries:
        if "_MASK" not in header:
            final_entries.append((header, sequence, structure))
        else:
            print(f"  Skipping MASK sequence from old file: {header}")

    # 2. new_mask_source_path의 모든 서열들을 추가
    #    (이 파일의 서열들은 모두 MASK 서열로 간주ตามคำอธิบาย)
    print(f"  Adding all sequences from new MASK source file ({len(new_mask_entries)} entries)")
    for header, sequence, structure in new_mask_entries:
        # 만약 new_mask_source_path 파일에도 _MASK가 아닌 것이 있을 가능성을 대비한다면, 여기서도 필터링 가능
        # if "_MASK" in header:
        #    final_entries.append((header, sequence, structure))
        # else:
        #    print(f"Warning: Non-MASK sequence found in new_mask_source_path and skipped: {header}")
        final_entries.append((header, sequence, structure))


    # 3. 최종 결과 파일 작성
    try:
        with open(output_fasta_path, 'w') as outfile:
            for header, sequence, structure in final_entries:
                outfile.write(f"{header}\n")
                outfile.write(f"{sequence}\n")
                outfile.write(f"{structure}\n")
        print(f"\nProcessing complete. Output written to: {output_fasta_path}")
        print(f"Total entries in output: {len(final_entries)}")
    except Exception as e:
        print(f"Error writing output file {output_fasta_path}: {e}")

if __name__ == "__main__":
    main()