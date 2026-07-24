#!/usr/bin/env python3
"""
pangenome_dispensable_windows.py

Computes dispensable node sharing between all paths in a pangenome variation
graph (GFA format) in sliding windows of X kb.

Dispensable nodes = nodes present in more than one but not all genomes.

Two modes:
  Default   : counts shared dispensable bp regardless of node order
  --colinear: finds co-linear blocks globally across the full path,
              then bins the co-linear bp into windows.
              Uses sparse DP - only processes focal positions that
              match each comp node, skipping the vast majority of
              work for typical pangenome dispensable node distributions.

Usage:
  python pangenome_dispensable_windows.py -i graph.gfa -w 100000 -o out.csv
  python pangenome_dispensable_windows.py -i graph.gfa -w 100000 -o out.csv --colinear -t 8
"""

import argparse
import os
from collections import defaultdict
from multiprocessing import current_process

import numpy as np
import pandas as pd


# =============================================================================
# GFA PARSING
# =============================================================================

def process_gfa_stream(file_path):
    print("Parsing GFA - pass 1: node lengths...")
    raw_lengths = {}
    with open(file_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip().split('\t')
            if parts[0] != 'S':
                continue
            raw_lengths[parts[1]] = len(parts[2])

    node_index   = {name: i for i, name in enumerate(raw_lengths)}
    node_lengths = np.zeros(len(node_index), dtype=np.int32)
    for name, length in raw_lengths.items():
        node_lengths[node_index[name]] = length
    del raw_lengths

    print("Parsing GFA - pass 2: paths...")
    paths = {}
    with open(file_path) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip().split('\t')
            if parts[0] != 'P':
                continue
            path_name  = parts[1]
            node_names = [n.rstrip('+-') for n in parts[2].split(',')]
            paths[path_name] = np.array(
                [node_index[n] for n in node_names], dtype=np.int32
            )

    print(f"  Nodes  : {len(node_index)}")
    print(f"  Paths  : {len(paths)}")
    print(f"  node_lengths array : {node_lengths.nbytes / 1e6:.1f} MB")
    print(f"  path arrays total  : "
          f"{sum(a.nbytes for a in paths.values()) / 1e6:.1f} MB")
    return node_lengths, node_index, paths


# =============================================================================
# NODE CLASSIFICATION
# =============================================================================

def classify_nodes(paths, n_nodes):
    print("Classifying nodes...")
    genome_name_to_id = {}
    genome_map        = {}
    for path_name in paths:
        genome_name = path_name.split('#')[0]
        if genome_name not in genome_name_to_id:
            genome_name_to_id[genome_name] = len(genome_name_to_id)
        genome_map[path_name] = genome_name_to_id[genome_name]

    n_genomes = len(genome_name_to_id)
    print(f"  Genomes detected : {n_genomes}")
    print(f"  Genome names     : {list(genome_name_to_id.keys())}")

    genome_node_sets = [set() for _ in range(n_genomes)]
    for path_name, node_ids in paths.items():
        gid = genome_map[path_name]
        genome_node_sets[gid].update(node_ids.tolist())

    node_genome_count = np.zeros(n_nodes, dtype=np.int16)
    for gid in range(n_genomes):
        ids = np.array(list(genome_node_sets[gid]), dtype=np.int32)
        node_genome_count[ids] += 1
        genome_node_sets[gid] = None

    dispensable = (node_genome_count > 1) & (node_genome_count < n_genomes)
    private     = node_genome_count == 1
    core        = node_genome_count == n_genomes

    print(f"  Core nodes        : {core.sum()}")
    print(f"  Dispensable nodes : {dispensable.sum()}")
    print(f"  Private nodes     : {private.sum()}")
    return dispensable


# =============================================================================
# WINDOW ACCUMULATION
# =============================================================================

def accumulate_overlaps(starts, ends, n_windows, window_size):
    counts = np.zeros(n_windows, dtype=np.int64)
    if len(starts) == 0:
        return counts

    first_win = starts // window_size
    last_win  = (ends - 1) // window_size

    for i in range(len(starts)):
        fw = int(first_win[i])
        lw = int(last_win[i])
        if fw == lw:
            counts[fw] += ends[i] - starts[i]
        else:
            counts[fw] += (fw + 1) * window_size - starts[i]
            if lw - fw > 1:
                counts[fw + 1: lw] += window_size
            counts[lw] += ends[i] - lw * window_size
    return counts


# =============================================================================
# POSITION COMPUTATION
# =============================================================================

def compute_disp_positions(nodes, node_lengths, disp_mask, chunk_size=500_000):
    """
    Compute genomic positions for dispensable nodes only.
    Gaps between consecutive dispensable nodes correctly reflect
    intervening core or private sequence.
    """
    n_nodes      = len(nodes)
    d_node_list  = []
    d_start_list = []
    d_end_list   = []
    pos          = np.int64(0)

    for chunk_start in range(0, n_nodes, chunk_size):
        chunk_end    = min(chunk_start + chunk_size, n_nodes)
        chunk_idx    = slice(chunk_start, chunk_end)
        chunk_nodes  = nodes[chunk_idx]
        chunk_lens   = node_lengths[chunk_nodes].astype(np.int64)
        chunk_disp   = disp_mask[chunk_idx]
        chunk_ends   = np.cumsum(chunk_lens) + pos
        chunk_starts = chunk_ends - chunk_lens

        if chunk_disp.any():
            d_node_list.append(chunk_nodes[chunk_disp])
            d_start_list.append(chunk_starts[chunk_disp])
            d_end_list.append(chunk_ends[chunk_disp])
        pos = chunk_ends[-1]

    if not d_node_list:
        return (
            np.array([], dtype=np.int32),
            np.array([], dtype=np.int64),
            np.array([], dtype=np.int64),
            int(pos),
        )
    return (
        np.concatenate(d_node_list),
        np.concatenate(d_start_list),
        np.concatenate(d_end_list),
        int(pos),
    )


# =============================================================================
# BUILD ORDERED DISPENSABLE SEQUENCES
# =============================================================================

def build_ordered_disp_sequences(paths, node_lengths, dispensable,
                                  chunk_size=500_000):
    print("Building ordered dispensable sequences...")
    disp_ordered = {}
    for path_name, nodes in paths.items():
        node_ids, starts, ends, path_len = compute_disp_positions(
            nodes, node_lengths, dispensable[nodes], chunk_size
        )
        disp_ordered[path_name] = {
            'node_ids' : node_ids,
            'starts'   : starts,
            'ends'     : ends,
            'path_len' : path_len,
        }
    total_mb = sum(
        d['node_ids'].nbytes + d['starts'].nbytes + d['ends'].nbytes
        for d in disp_ordered.values()
    ) / 1e6
    print(f"  disp_ordered total : {total_mb:.1f} MB")
    return disp_ordered


# =============================================================================
# CO-LINEAR BLOCK DETECTION - SPARSE DP
#
# Key insight: for a typical pangenome, each dispensable node appears
# in only a few paths out of all paths. This means for most comp positions
# ci, very few focal positions will match. The dense DP wastes time
# checking all n_focal positions even when only a handful match.
#
# Sparse approach:
#   1. Build focal_index: node_id -> list of positions in focal
#      Built once per comparison, O(n_focal)
#   2. For each comp position ci:
#      - look up focal positions matching comp_ids[ci]
#      - only process those positions (typically ~10, not 650K)
#      - skip entirely if no matches
#   3. Track active runs in a dict (prev_runs: fi -> run_length)
#      Only positions currently in a run are stored
#
# Complexity:
#   Dense:  O(n_focal * n_comp) numpy ops
#   Sparse: O(n_comp * avg_matches) Python ops
#           avg_matches = n_focal / n_unique_nodes * avg_copies_per_node
#           for typical pangenome: avg_matches ~10-100
#           speedup: 650K / 10 = ~65000x
#
# Memory:
#   Dense:  O(n_focal) arrays always fully allocated
#   Sparse: prev_runs dict only stores active run positions
#           typically much smaller than n_focal
# =============================================================================

def find_colinear_blocks_with_repeats(focal_ids, focal_starts, focal_ends,
                                       comp_ids,  comp_starts,  comp_ends,
                                       min_block_size=2):
    """
    Sparse co-linear block detection using inverted index.

    For each comp position, only processes focal positions that actually
    match the comp node. Skips columns with no matches entirely.

    Parameters
    ----------
    focal_ids    : int32 array
    focal_starts : int64 array
    focal_ends   : int64 array
    comp_ids     : int32 array
    comp_starts  : int64 array
    comp_ends    : int64 array
    min_block_size : int

    Returns
    -------
    focal_mask : bool array over focal_ids
    """
    n_focal    = len(focal_ids)
    n_comp     = len(comp_ids)
    focal_mask = np.zeros(n_focal, dtype=bool)

    if n_focal == 0 or n_comp == 0:
        return focal_mask

    # Precompute focal adjacency once
    # focal_adj[fi] = True if no gap between focal node fi-1 and fi
    focal_gaps     = np.empty(n_focal, dtype=np.int64)
    focal_gaps[0]  = -1
    focal_gaps[1:] = focal_starts[1:] - focal_ends[:-1]
    focal_adj      = (focal_gaps == 0)

    # Precompute comp adjacency once
    comp_gaps      = np.empty(n_comp, dtype=np.int64)
    comp_gaps[0]   = -1
    comp_gaps[1:]  = comp_starts[1:] - comp_ends[:-1]

    # Build inverted index: node_id -> sorted list of focal positions
    # This lets us jump directly to matching focal positions for each ci
    focal_index = defaultdict(list)
    for fi in range(n_focal):
        focal_index[int(focal_ids[fi])].append(fi)

    # Best block tracking
    # best_len[fi]   = length of longest completed block ending at fi
    # best_start[fi] = start index of that block
    best_len   = np.zeros(n_focal, dtype=np.int32)
    best_start = np.arange(n_focal, dtype=np.int32)

    # Sparse active runs: fi -> run_length
    # Only positions currently in an active run are stored
    # Much smaller than a full n_focal array when matches are sparse
    prev_runs = {}

    for ci in range(n_comp):
        nid         = int(comp_ids[ci])
        comp_adj_ci = bool(comp_gaps[ci] == 0)

        # Get focal positions matching this comp node
        matching_fi = focal_index.get(nid)

        if not matching_fi:
            # No focal positions match this comp node
            # All active runs that cannot continue must be completed
            if prev_runs:
                for fi, run_len in prev_runs.items():
                    if run_len >= min_block_size:
                        fi_start = fi - run_len + 1
                        if run_len > int(best_len[fi]):
                            best_len[fi]   = run_len
                            best_start[fi] = fi_start
                prev_runs = {}
            continue

        # Build new active runs for this column
        curr_runs = {}

        for fi in matching_fi:
            prev_fi  = fi - 1
            run_len  = prev_runs.get(prev_fi, 0)

            if run_len > 0 and bool(focal_adj[fi]) and comp_adj_ci:
                # Extend existing run
                curr_runs[fi] = run_len + 1
            else:
                # Start new run of length 1
                curr_runs[fi] = 1

        # Find runs in prev_runs that were NOT continued
        # A run at fi continues if curr_runs contains fi+1
        # with length == prev run length + 1
        for fi, run_len in prev_runs.items():
            next_fi = fi + 1
            if curr_runs.get(next_fi, 0) != run_len + 1:
                # Run ended at fi
                if run_len >= min_block_size:
                    fi_start = fi - run_len + 1
                    if run_len > int(best_len[fi]):
                        best_len[fi]   = run_len
                        best_start[fi] = fi_start

        prev_runs = curr_runs

    # Handle runs still open after last column
    for fi, run_len in prev_runs.items():
        if run_len >= min_block_size:
            fi_start = fi - run_len + 1
            if run_len > int(best_len[fi]):
                best_len[fi]   = run_len
                best_start[fi] = fi_start

    # Select non-overlapping blocks greedily by length descending
    valid = best_len >= min_block_size
    if not valid.any():
        return focal_mask

    valid_fi   = np.where(valid)[0]
    sort_order = np.argsort(-best_len[valid_fi])
    sorted_fi  = valid_fi[sort_order]

    focal_used = np.zeros(n_focal, dtype=bool)
    for fi in sorted_fi:
        length      = int(best_len[fi])
        fi_start    = int(best_start[fi])
        focal_range = slice(fi_start, fi_start + length)
        if focal_used[focal_range].any():
            continue
        focal_mask[focal_range] = True
        focal_used[focal_range] = True

    return focal_mask


# =============================================================================
# PROCESS-LEVEL GLOBALS
# =============================================================================

GLOBAL_NODE_LENGTHS = None
GLOBAL_DISPENSABLE  = None
GLOBAL_PATH_NODES   = None
GLOBAL_DISP_SETS    = None
GLOBAL_DISP_ORDERED = None
GLOBAL_WINDOW       = None
GLOBAL_OUTPUT       = None


# =============================================================================
# WORKER: UNORDERED MODE
# =============================================================================

def process_focal_unordered(focal):
    node_lengths = GLOBAL_NODE_LENGTHS
    dispensable  = GLOBAL_DISPENSABLE
    path_nodes   = GLOBAL_PATH_NODES
    disp_sets    = GLOBAL_DISP_SETS
    window_size  = GLOBAL_WINDOW
    output_base  = GLOBAL_OUTPUT

    f_genome = focal.split('#')[0]
    nodes    = path_nodes[focal]

    disp_mask = dispensable[nodes]
    if not disp_mask.any():
        return

    d_node_ids, d_starts, d_ends, focal_length = compute_disp_positions(
        nodes, node_lengths, disp_mask
    )

    n_windows = (focal_length - 1) // window_size + 1
    focal_bp  = accumulate_overlaps(d_starts, d_ends, n_windows, window_size)

    window_lengths     = np.full(n_windows, window_size, dtype=np.int64)
    window_lengths[-1] = focal_length - (n_windows - 1) * window_size

    pid          = current_process().pid
    outfile      = f"{output_base}.part_{pid}.csv"
    write_header = not os.path.exists(outfile)

    for comp, comp_disp_ids in disp_sets.items():
        if comp == focal:
            continue
        if comp.split('#')[0] == f_genome:
            continue

        shared_mask = np.isin(d_node_ids, comp_disp_ids, assume_unique=False)
        if not shared_mask.any():
            continue

        shared_bp = accumulate_overlaps(
            d_starts[shared_mask],
            d_ends[shared_mask],
            n_windows,
            window_size,
        )

        active = focal_bp > 0
        if not active.any():
            continue

        w_idx    = np.where(active)[0]
        pct      = (shared_bp[w_idx] / window_lengths[w_idx]) * 100.0
        b_starts = w_idx * window_size + 1
        b_ends   = np.minimum((w_idx + 1) * window_size, focal_length)

        df = pd.DataFrame({
            "Chromosome_Window_Start" : b_starts,
            "Chromosome_Window_End"   : b_ends,
            "Focal_Path"              : focal,
            "Compared_Path"           : comp,
            "Shared_pct"              : pct,
            "Shared_bp"               : shared_bp[w_idx],
            "Focal_disp_bp"           : focal_bp[w_idx],
            "Window_length_bp"        : window_lengths[w_idx],
        })
        df.to_csv(outfile, mode='a', header=write_header, index=False)
        write_header = False


# =============================================================================
# WORKER: CO-LINEAR MODE
# One task per (focal, comp) pair - memory released after each comparison
# =============================================================================

def process_focal_comp_pair(args):
    """
    Process a single focal vs comp comparison.
    One task = one pair = memory fully released after each call.
    """
    focal, comp  = args
    disp_ordered = GLOBAL_DISP_ORDERED
    window_size  = GLOBAL_WINDOW
    output_base  = GLOBAL_OUTPUT

    focal_data = disp_ordered[focal]
    comp_data  = disp_ordered[comp]

    if len(focal_data['node_ids']) == 0:
        return
    if len(comp_data['node_ids']) == 0:
        return

    focal_ids    = focal_data['node_ids']
    focal_starts = focal_data['starts']
    focal_ends   = focal_data['ends']
    focal_length = focal_data['path_len']

    n_windows = (focal_length - 1) // window_size + 1

    focal_bp = accumulate_overlaps(
        focal_starts, focal_ends, n_windows, window_size
    )

    window_lengths     = np.full(n_windows, window_size, dtype=np.int64)
    window_lengths[-1] = focal_length - (n_windows - 1) * window_size

    colinear_mask = find_colinear_blocks_with_repeats(
        focal_ids,
        focal_starts,
        focal_ends,
        comp_data['node_ids'],
        comp_data['starts'],
        comp_data['ends'],
    )

    if not colinear_mask.any():
        return

    shared_bp = accumulate_overlaps(
        focal_starts[colinear_mask],
        focal_ends[colinear_mask],
        n_windows,
        window_size,
    )

    active = focal_bp > 0
    if not active.any():
        return

    w_idx    = np.where(active)[0]
    pct      = (shared_bp[w_idx] / window_lengths[w_idx]) * 100.0
    b_starts = w_idx * window_size + 1
    b_ends   = np.minimum((w_idx + 1) * window_size, focal_length)

    pid          = current_process().pid
    outfile      = f"{output_base}.part_{pid}.csv"
    write_header = not os.path.exists(outfile)

    df = pd.DataFrame({
        "Chromosome_Window_Start" : b_starts,
        "Chromosome_Window_End"   : b_ends,
        "Focal_Path"              : focal,
        "Compared_Path"           : comp,
        "Shared_pct"              : pct,
        "Shared_bp"               : shared_bp[w_idx],
        "Focal_disp_bp"           : focal_bp[w_idx],
        "Window_length_bp"        : window_lengths[w_idx],
    })
    df.to_csv(outfile, mode='a', header=write_header, index=False)


# =============================================================================
# MAIN CALCULATION DISPATCHER
# =============================================================================

def calculate(paths, node_lengths, dispensable, window_size, output,
              n_workers, colinear=False):
    global GLOBAL_NODE_LENGTHS, GLOBAL_DISPENSABLE, GLOBAL_PATH_NODES
    global GLOBAL_DISP_SETS, GLOBAL_DISP_ORDERED, GLOBAL_WINDOW, GLOBAL_OUTPUT

    GLOBAL_NODE_LENGTHS = node_lengths
    GLOBAL_DISPENSABLE  = dispensable
    GLOBAL_PATH_NODES   = paths
    GLOBAL_WINDOW       = window_size
    GLOBAL_OUTPUT       = output

    import multiprocessing
    ctx = multiprocessing.get_context('fork')

    if colinear:
        disp_ordered        = build_ordered_disp_sequences(
            paths, node_lengths, dispensable
        )
        GLOBAL_DISP_ORDERED = disp_ordered

        path_list = list(paths.keys())
        pairs = [
            (focal, comp)
            for focal in path_list
            for comp  in path_list
            if focal != comp
            and focal.split('#')[0] != comp.split('#')[0]
        ]
        print(f"Mode: co-linear  |  workers: {n_workers}  "
              f"|  pairs: {len(pairs)}")

        pool = ctx.Pool(processes=n_workers)
        try:
            for i, _ in enumerate(
                pool.imap_unordered(process_focal_comp_pair, pairs), 1
            ):
                print(f"  Finished {i}/{len(pairs)}", flush=True)
        finally:
            pool.close()
            pool.join()

    else:
        print("Building dispensable ID sets (unordered mode)...")
        disp_sets     = {}
        total_disp_mb = 0
        for p, nodes in paths.items():
            d              = np.unique(nodes[dispensable[nodes]]).astype(np.int32)
            disp_sets[p]   = d
            total_disp_mb += d.nbytes
        print(f"  disp_sets total : {total_disp_mb / 1e6:.1f} MB")
        GLOBAL_DISP_SETS = disp_sets

        path_list = list(paths.keys())
        print(f"Mode: unordered  |  workers: {n_workers}")

        pool = ctx.Pool(processes=n_workers)
        try:
            for i, _ in enumerate(
                pool.imap_unordered(process_focal_unordered, path_list), 1
            ):
                print(f"  Finished {i}/{len(path_list)}", flush=True)
        finally:
            pool.close()
            pool.join()


# =============================================================================
# ENTRY POINT
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Compute dispensable node sharing between pangenome paths "
            "in sliding windows."
        )
    )
    parser.add_argument("-i", "--input",   required=True, metavar="GFA")
    parser.add_argument("-w", "--window",  required=True, type=int, metavar="BP")
    parser.add_argument("-o", "--output",  default="output.csv", metavar="CSV")
    parser.add_argument("-t", "--threads", default=4, type=int, metavar="N")
    parser.add_argument(
        "--colinear",
        action="store_true",
        help=(
            "Find co-linear blocks using sparse DP. One task per "
            "(focal,comp) pair. Memory released between comparisons."
        ),
    )
    parser.add_argument(
        "--min-block-size",
        default=2,
        type=int,
        metavar="N",
        help="Minimum co-linear block size in nodes (default: 2).",
    )

    args = parser.parse_args()

    node_lengths, node_index, paths = process_gfa_stream(args.input)
    dispensable = classify_nodes(paths, len(node_lengths))
    calculate(
        paths, node_lengths, dispensable,
        args.window, args.output, args.threads,
        colinear=args.colinear,
    )

    print("\nDone.")
    print("Merge part files with:")
    print(
        f"  awk 'NR==1 || !/^Chromosome_Window_Start/' "
        f"{args.output}.part_*.csv > {args.output}"
    )
    print("Then clean up with:")
    print(f"  rm {args.output}.part_*.csv")


if __name__ == "__main__":
    main()
