#!/usr/bin/env python3
"""
pangenome_dispensable_windows_debug.py

Debug version that prints co-linear blocks found between paths.

Usage:
  python diagnostic.py -i graph.gfa
  python diagnostic.py -i graph.gfa --focal DroSim#1#2L
  python diagnostic.py -i graph.gfa --focal DroSim#1#2L --min-block-size 3
"""

import argparse
import numpy as np


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
    index_node   = {i: name for name, i in node_index.items()}
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
    return node_lengths, node_index, index_node, paths


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
    core        = node_genome_count == n_genomes
    private     = node_genome_count == 1

    print(f"  Core nodes        : {core.sum()}")
    print(f"  Dispensable nodes : {dispensable.sum()}")
    print(f"  Private nodes     : {private.sum()}")
    return dispensable


# =============================================================================
# POSITION COMPUTATION
# =============================================================================

def compute_disp_positions(nodes, node_lengths, disp_mask, chunk_size=500_000):
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
    return disp_ordered


# =============================================================================
# CO-LINEAR BLOCK DETECTION - VECTORISED
# =============================================================================

def find_colinear_blocks_with_repeats(focal_ids, focal_starts, focal_ends,
                                       comp_ids,  comp_starts,  comp_ends,
                                       min_block_size=2):
    n_focal    = len(focal_ids)
    n_comp     = len(comp_ids)
    focal_mask = np.zeros(n_focal, dtype=bool)

    if n_focal == 0 or n_comp == 0:
        return focal_mask

    focal_gaps     = np.empty(n_focal, dtype=np.int64)
    focal_gaps[0]  = -1
    focal_gaps[1:] = focal_starts[1:] - focal_ends[:-1]
    focal_adj      = (focal_gaps == 0)

    comp_gaps      = np.empty(n_comp, dtype=np.int64)
    comp_gaps[0]   = -1
    comp_gaps[1:]  = comp_starts[1:] - comp_ends[:-1]

    prev_col     = np.zeros(n_focal, dtype=np.int32)
    curr_col     = np.zeros(n_focal, dtype=np.int32)
    prev_shifted = np.zeros(n_focal, dtype=np.int32)
    all_blocks   = []

    for ci in range(n_comp):
        match = (focal_ids == comp_ids[ci])

        prev_shifted[0]  = 0
        prev_shifted[1:] = prev_col[:-1]

        comp_adj_ci = (comp_gaps[ci] == 0)
        can_extend  = (prev_shifted > 0) & focal_adj & comp_adj_ci

        extend  = match & can_extend
        new_run = match & ~extend

        curr_col[:]       = 0
        curr_col[extend]  = prev_shifted[extend] + 1
        curr_col[new_run] = 1

        run_ends = prev_col >= min_block_size
        if run_ends.any():
            continues      = np.zeros(n_focal, dtype=bool)
            continues[:-1] = (curr_col[1:] == prev_col[:-1] + 1)
            completed      = run_ends & ~continues
            if completed.any():
                for fi in np.where(completed)[0]:
                    run_len  = int(prev_col[fi])
                    fi_start = fi - run_len + 1
                    all_blocks.append((run_len, fi_start))

        prev_col, curr_col = curr_col, prev_col

    run_ends = prev_col >= min_block_size
    if run_ends.any():
        for fi in np.where(run_ends)[0]:
            run_len  = int(prev_col[fi])
            fi_start = fi - run_len + 1
            all_blocks.append((run_len, fi_start))

    if not all_blocks:
        return focal_mask

    all_blocks.sort(key=lambda x: -x[0])
    focal_used = np.zeros(n_focal, dtype=bool)

    for length, fi_start in all_blocks:
        focal_range = slice(fi_start, fi_start + length)
        if focal_used[focal_range].any():
            continue
        focal_mask[focal_range] = True
        focal_used[focal_range] = True

    return focal_mask


# =============================================================================
# PRINT CO-LINEAR BLOCKS
# =============================================================================

def print_colinear_blocks(focal_ids, focal_starts, focal_ends,
                           focal_mask, index_node):
    """
    Print all co-linear blocks found in focal_mask as contiguous runs.
    """
    if not focal_mask.any():
        print("    No co-linear blocks found")
        return

    # Find contiguous runs of True in focal_mask
    in_block    = False
    block_start = 0
    blocks      = []

    for i in range(len(focal_mask)):
        if focal_mask[i] and not in_block:
            block_start = i
            in_block    = True
        elif not focal_mask[i] and in_block:
            blocks.append((block_start, i))
            in_block = False
    if in_block:
        blocks.append((block_start, len(focal_mask)))

    total_bp = 0
    for b_idx, (b_start, b_end) in enumerate(blocks):

        node_names = [index_node[int(nid)] for nid in focal_ids[b_start:b_end]]
        start_pos  = int(focal_starts[b_start])
        end_pos    = int(focal_ends[b_end - 1])
        block_bp   = int(np.sum(
            focal_ends[b_start:b_end] - focal_starts[b_start:b_end]
        ))
        total_bp  += block_bp

        print(f"    Block {b_idx + 1}:")
        print(f"      n_nodes  : {b_end - b_start}")
        print(f"      bp       : {block_bp}")
        print(f"      position : [{start_pos} - {end_pos})")
        print(f"      nodes    : {', '.join(node_names)}")

    print(f"    ----")
    print(f"    Total blocks       : {len(blocks)}")
    print(f"    Total co-linear bp : {total_bp}")


# =============================================================================
# MAIN DEBUG FUNCTION
# =============================================================================

def run_debug(file_path, focal_filter=None, min_block_size=2):

    node_lengths, node_index, index_node, paths = process_gfa_stream(file_path)
    dispensable  = classify_nodes(paths, len(node_lengths))
    disp_ordered = build_ordered_disp_sequences(paths, node_lengths,
                                                dispensable)

    print()
    print("=" * 70)
    print("CO-LINEAR BLOCK ANALYSIS")
    print("=" * 70)

    for focal_name, focal_data in disp_ordered.items():

        # Filter to requested focal path if specified
        if focal_filter and focal_filter not in focal_name:
            continue

        f_genome     = focal_name.split('#')[0]
        focal_ids    = focal_data['node_ids']
        focal_starts = focal_data['starts']
        focal_ends   = focal_data['ends']
        focal_length = focal_data['path_len']

        print()
        print(f"FOCAL: {focal_name}")
        print(f"  Path length       : {focal_length} bp")
        print(f"  Dispensable nodes : {len(focal_ids)}")

        if len(focal_ids) == 0:
            print("  No dispensable nodes, skipping")
            continue

        # Print focal dispensable sequence
        print("  Dispensable sequence:")
        for i in range(len(focal_ids)):
            node_name = index_node[int(focal_ids[i])]
            start     = int(focal_starts[i])
            end       = int(focal_ends[i])
            bp        = end - start

            if i > 0:
                gap = int(focal_starts[i]) - int(focal_ends[i - 1])
                gap_str = f"  gap={gap}" if gap > 0 else ""
            else:
                gap_str = ""

            print(f"    [{i:5d}]  node={node_name:>15s}"
                  f"  pos=[{start:10d}, {end:10d})"
                  f"  bp={bp:8d}"
                  f"{gap_str}")

        print()

        # Compare against every other path
        for comp_name, comp_data in disp_ordered.items():
            if comp_name == focal_name:
                continue
            if comp_name.split('#')[0] == f_genome:
                continue
            if len(comp_data['node_ids']) == 0:
                continue

            print(f"  vs COMP: {comp_name}")
            print(f"    Comp dispensable nodes : {len(comp_data['node_ids'])}")

            colinear_mask = find_colinear_blocks_with_repeats(
                focal_ids,
                focal_starts,
                focal_ends,
                comp_data['node_ids'],
                comp_data['starts'],
                comp_data['ends'],
                min_block_size=min_block_size,
            )

            print_colinear_blocks(
                focal_ids,
                focal_starts,
                focal_ends,
                colinear_mask,
                index_node,
            )
            print()

        print("-" * 70)


# =============================================================================
# ENTRY POINT
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Print co-linear dispensable blocks between paths."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        metavar="GFA",
        help="Input GFA file.",
    )
    parser.add_argument(
        "--focal",
        default=None,
        metavar="NAME",
        help=(
            "Only process paths containing this string as focal. "
            "e.g. --focal DroSim#1#2L processes all paths matching "
            "that string. If not specified all paths are processed."
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
    run_debug(args.input, args.focal, args.min_block_size)


if __name__ == "__main__":
    main()
