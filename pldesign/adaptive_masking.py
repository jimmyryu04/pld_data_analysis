#!/usr/bin/env python3
#
# Copyright (c) 2024-2025 Seoul National University
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#

import numpy as np
from collections import deque

def find_antidiagonal_groups(matrix, max_gap_size, min_group_size=5):
    if not isinstance(matrix, np.ndarray) or matrix.ndim != 2:
        raise ValueError("Input 'grid' must be a 2D NumPy array.")
    if not isinstance(max_gap_size, int) or max_gap_size < 0:
        raise ValueError("'max_gap_size' must be a non-negative integer.")
    if not isinstance(min_group_size, int) or min_group_size < 1:
        raise ValueError("'min_group_size' must be a positive integer.")

    rows, cols = matrix.shape
    visited = set() # Keep track of cells already assigned to a group
    groups = []
    shift_allowance = max_gap_size

    # Define anti-diagonal directions
    directions = [(-1, 1), (1, -1)] # (dr, dc) for top-right/bottom-left

    # Helper to check if a cell is within grid bounds
    def is_within_bounds(r, c):
        return 0 <= r < rows and 0 <= c < cols

    # Helper to check if a cell is active (True or non-zero)
    def is_active(r, c):
        # Check bounds first to avoid index errors
        if is_within_bounds(r, c):
            # Handles both boolean and integer grids
            return matrix[r, c] != 0
        return False

    # Iterate through each cell in the grid
    for r_start in range(rows):
        for c_start in range(cols):
            # If the cell is active and hasn't been visited as part of another group
            if is_active(r_start, c_start) and (r_start, c_start) not in visited:
                # Start a Breadth-First Search (BFS) to find the connected group
                current_group = []
                q = deque([(r_start, c_start)])
                visited.add((r_start, c_start)) # Mark starting cell as visited
                current_group.append((r_start, c_start))

                while q:
                    r, c = q.popleft()

                    # --- Find Neighbors (Direct and Gap-Bridged) ---

                    potential_neighbors = []

                    # 1. Check Direct Anti-diagonal Neighbors
                    for dr, dc in directions:
                        nr, nc = r + dr, c + dc
                        if is_active(nr, nc):
                            potential_neighbors.append((nr, nc))

                    # 2. Check Gap Bridging (if allowed)
                    if max_gap_size > 0:
                        for dr, dc in directions:
                            for g in range(1, max_gap_size + 1):
                                # Check if path for gap 'g' is clear of active cells
                                path_clear = True
                                for i in range(1, g + 1):
                                    ir, ic = r + i * dr, c + i * dc
                                    # If path goes out of bounds or hits an active cell, it's blocked
                                    if not is_within_bounds(ir, ic) or is_active(ir, ic):
                                        path_clear = False
                                        break # Stop checking this path

                                if not path_clear:
                                    # Cannot bridge further in this direction if path is blocked
                                    break # Move to next direction or stop checking larger gaps

                                # Path of size 'g' is clear, check the target area
                                # Calculate ideal target position P
                                pr, pc = r + (g + 1) * dr, c + (g + 1) * dc

                                # Check if ideal target P is itself out of bounds; if so, no need to search window
                                if not is_within_bounds(pr, pc):
                                     # Although P is out of bounds, the search window might still overlap
                                     # with the grid. Continue to search window calculation.
                                     pass # Let the window search handle bounds.

                                # Define the search window centered at P
                                max_r_search = pr + shift_allowance
                                min_c_search = pc - shift_allowance

                                # Search within the window for *any* active cell
                                for sr in range(pr, max_r_search + 1):
                                    for sc in range(min_c_search, pc + 1):
                                        # Check if the cell in the search window is active
                                        if is_active(sr, sc):
                                             # Found an active cell B in the search window
                                             potential_neighbors.append((sr, sc))

                                # Optimization: If we found a connection across gap 'g',
                                # we don't need to check larger gaps (g+1, ...) in the same direction
                                # from the *current* cell (r, c), because the path check for g+1
                                # would fail at the g+1 step anyway if P was active. If P wasn't active
                                # but some other cell (sr, sc) in the window *was* active, the BFS
                                # will eventually explore from (sr, sc).
                                # However, the rule connects A to *any* B found in the window for that gap size g.
                                # The current logic correctly finds all potential neighbors via gaps.
                                # Let's stick with the path_clear check as the primary way to stop increasing g.


                    # --- Process Neighbors ---
                    for nr, nc in potential_neighbors:
                        # If the neighbor is within bounds (redundant check?) and not visited
                        # is_active already checks bounds implicitly for finding neighbors
                        if (nr, nc) not in visited:
                            visited.add((nr, nc))
                            current_group.append((nr, nc))
                            q.append((nr, nc))

                # BFS for this group is complete, check size requirement
                if len(current_group) >= min_group_size:
                    groups.append(current_group)

    return groups

def dotbracket_to_matrix(dot_bracket):
    seq_len = len(dot_bracket)

    bp_matrix = np.zeros((seq_len, seq_len), dtype=np.int8)

    open_indices_stack = []

    for i, notation in enumerate(dot_bracket):
        if notation == '(':
            open_indices_stack.append(i)
        elif notation == ')':
            if not open_indices_stack:
                raise ValueError('Invalid dot-bracket string: Unmatched closing '
                                 f"bracket ')' at index {i}.")
            j = open_indices_stack.pop()

            bp_matrix[i, j] = 1
        elif notation == '.':
            pass
        else:
             raise ValueError(f"Invalid character '{notation}' in dot-bracket "
                              f"string at index {i}. Only '(', ')', '.' are allowed.")

    if open_indices_stack:
        unmatched_indices = ', '.join(map(str, open_indices_stack))
        raise ValueError("Invalid dot-bracket string: Unmatched opening bracket(s) "
                         f"'(' at indices: {unmatched_indices}.")

    return bp_matrix

def find_masking_regions(structure, max_gap_size=2, min_group_size=5, max_group_size=8):
    matrix = dotbracket_to_matrix(structure)
    groups = find_antidiagonal_groups(matrix, max_gap_size, min_group_size)

    masking_regions = []
    for group in groups:
        num_subregions = (len(group) + max_group_size - 1) // max_group_size
        edges = np.linspace(0, len(group), num_subregions + 1, endpoint=True, dtype=int)
        for left, right in zip(edges[:-1], edges[1:]):
            masking_regions.append(group[left:right])

    return sorted([((r[-1][1], r[0][1]+1), (r[0][0], r[-1][0]+1)) for r in masking_regions])

if __name__ == "__main__":
    import pandas as pd

    results = pd.read_csv('eval-gamma-lambda/result-free.csv')

    min_mask_regions = 0.1
    max_mask_regions = 0.3
    num_mask_groups = 10

    regions = find_masking_regions(results.iloc[0]['structure'], max_gap_size=2, min_group_size=5, max_group_size=8)
    regionindices = list(range(len(regions)))

    masks = []

    for i in range(num_mask_groups):
        num_masks = int(np.random.uniform(min_mask_regions, max_mask_regions) * len(regions))
        selected_masks = np.random.choice(regionindices, num_masks, replace=False)
        selected_masks = sorted(selected_masks)
        selected_masks = [regions[i] for i in selected_masks]
        masks.append(selected_masks)

    print(masks)
