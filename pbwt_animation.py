#!/usr/bin/env python3
"""
PBWT Animation - Visualizing Richard Durbin's Positional Burrows-Wheeler Transform

This script creates an animation showing how haplotype sequences are progressively
sorted by their reversed prefixes at each position, as shown in Figure 1 of:
Durbin R. "Efficient haplotype matching and storage using the positional
Burrows-Wheeler transform (PBWT)" Bioinformatics 2014.

The animation shows:
- Haplotypes sorted by reversed prefix at each position k
- The allele values at position k highlighted
- How the sort order changes from position k to k+1
- Maximal matches between adjacent sequences underlined
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.colors import LinearSegmentedColormap
import argparse


def simulate_haplotypes(n_haplotypes=8, n_sites=12, seed=42):
    """
    Simulate a small haplotype panel with some structure.
    Creates correlated haplotypes to show interesting sorting behavior.
    """
    np.random.seed(seed)

    # Create base haplotypes with some structure
    haplotypes = np.zeros((n_haplotypes, n_sites), dtype=int)

    # Create a few "founder" patterns and derive others from them
    n_founders = 3
    founders = np.random.randint(0, 2, (n_founders, n_sites))

    for i in range(n_haplotypes):
        # Pick a founder to base this haplotype on
        founder_idx = i % n_founders
        haplotypes[i] = founders[founder_idx].copy()

        # Add some random mutations
        n_mutations = np.random.randint(1, n_sites // 3)
        mutation_sites = np.random.choice(n_sites, n_mutations, replace=False)
        haplotypes[i, mutation_sites] = 1 - haplotypes[i, mutation_sites]

    return haplotypes


def compute_pbwt(haplotypes):
    """
    Compute the PBWT (positional prefix arrays) for all positions.

    Returns:
        a_list: List of permutation arrays (sort order at each position)
        d_list: List of divergence arrays (position where match with previous starts)
    """
    M, N = haplotypes.shape  # M haplotypes, N sites

    a_list = []  # Sort order at each position
    d_list = []  # Divergence values

    # Initial order (identity permutation)
    a = np.arange(M)
    d = np.zeros(M, dtype=int)

    a_list.append(a.copy())
    d_list.append(d.copy())

    for k in range(N):
        # Get allele values at position k in current sorted order
        y = haplotypes[a, k]

        # Build new sorted order by stable partitioning:
        # First all sequences with y[i]=0, then all with y[i]=1
        # This is the key insight: sorting by prefix up to k+1 just requires
        # splitting by value at k while preserving relative order

        a_new = np.empty(M, dtype=int)
        d_new = np.empty(M, dtype=int)

        # Track where 0s and 1s go
        zeros_idx = 0
        ones_idx = np.sum(y == 0)

        p = k + 1  # Track divergence for zeros group
        q = k + 1  # Track divergence for ones group

        for i in range(M):
            if y[i] == 0:
                a_new[zeros_idx] = a[i]
                # Divergence: if previous in new order had different value, match starts at k
                d_new[zeros_idx] = p
                p = d[i]
                zeros_idx += 1
            else:
                a_new[ones_idx] = a[i]
                d_new[ones_idx] = q
                q = d[i]
                ones_idx += 1

        a = a_new
        d = d_new
        a_list.append(a.copy())
        d_list.append(d.copy())

    return a_list, d_list


def create_figure1_frame(haplotypes, a_list, d_list, k, ax, original_indices=None):
    """
    Create a single frame showing the PBWT state at position k.

    This mimics Figure 1 from Durbin's paper:
    - Shows haplotypes sorted by reversed prefix up to position k
    - Highlights the current position k
    - Shows match lengths with previous sequence
    """
    ax.clear()

    M, N = haplotypes.shape
    a_k = a_list[k]  # Sort order at position k
    d_k = d_list[k]  # Divergence values at position k

    # Colors for alleles
    color_0 = '#4A90D9'  # Blue for 0
    color_1 = '#E74C3C'  # Red for 1
    bg_color = '#F5F5F5'
    highlight_color = '#FFE066'
    match_color = '#90EE90'  # Light green for matching regions

    cell_width = 1.0
    cell_height = 0.8

    # Draw each haplotype in sorted order
    for row_idx, hap_idx in enumerate(a_k):
        y_pos = M - row_idx - 1  # Flip so row 0 is at top

        # Get the haplotype sequence
        hap = haplotypes[hap_idx]

        # Calculate match length with previous sequence in sorted order
        if row_idx > 0:
            prev_hap_idx = a_k[row_idx - 1]
            match_start = d_k[row_idx]
        else:
            match_start = k  # No previous sequence

        # Draw each site
        for site in range(N):
            x_pos = site * cell_width

            # Determine cell color and style
            allele = hap[site]

            # Background color
            if site == k:
                # Highlight current position
                bg = highlight_color
            elif site < k and site >= match_start and row_idx > 0:
                # Part of match with previous sequence
                bg = match_color
            else:
                bg = bg_color

            # Draw cell background
            rect = patches.Rectangle(
                (x_pos, y_pos * cell_height),
                cell_width * 0.95,
                cell_height * 0.9,
                facecolor=bg,
                edgecolor='gray',
                linewidth=0.5
            )
            ax.add_patch(rect)

            # Draw allele - always visible, but faded for future positions
            if site <= k:
                # Processed positions: bold colors
                allele_color = color_1 if allele == 1 else color_0
                fontweight = 'bold'
                alpha = 1.0
            else:
                # Future positions: faded colors
                allele_color = '#E8A0A0' if allele == 1 else '#A0C4E8'
                fontweight = 'normal'
                alpha = 0.6

            ax.text(
                x_pos + cell_width * 0.45,
                y_pos * cell_height + cell_height * 0.45,
                str(allele),
                ha='center', va='center',
                fontsize=12, fontweight=fontweight,
                color=allele_color, alpha=alpha
            )

        # Add row labels (original haplotype index)
        ax.text(
            -0.5, y_pos * cell_height + cell_height * 0.45,
            f'h{hap_idx}',
            ha='right', va='center',
            fontsize=10, fontweight='bold'
        )

        # Add divergence value on the right
        if row_idx > 0:
            ax.text(
                N * cell_width + 0.3,
                y_pos * cell_height + cell_height * 0.45,
                f'd={match_start}',
                ha='left', va='center',
                fontsize=8,
                color='gray'
            )

    # Add position labels at bottom
    for site in range(N):
        color = 'red' if site == k else 'black'
        weight = 'bold' if site == k else 'normal'
        ax.text(
            site * cell_width + cell_width * 0.45,
            -0.4,
            str(site),
            ha='center', va='center',
            fontsize=9,
            color=color,
            fontweight=weight
        )

    # Set axis properties
    ax.set_xlim(-1, N * cell_width + 1.5)
    ax.set_ylim(-0.8, M * cell_height + 0.3)
    ax.set_aspect('equal')
    ax.axis('off')

    # Title
    ax.set_title(
        f'PBWT Visualization: Position k = {k}\n'
        f'Haplotypes sorted by reversed prefix [0..{k-1}]' if k > 0 else
        f'PBWT Visualization: Position k = {k}\n'
        f'Initial order (no prefix yet)',
        fontsize=12, fontweight='bold', pad=10
    )

    # Add legend
    legend_y = M * cell_height + 0.1
    ax.add_patch(patches.Rectangle((0, legend_y), 0.3, 0.2, facecolor=highlight_color, edgecolor='gray'))
    ax.text(0.5, legend_y + 0.1, 'Current position', fontsize=8, va='center')

    ax.add_patch(patches.Rectangle((4, legend_y), 0.3, 0.2, facecolor=match_color, edgecolor='gray'))
    ax.text(4.5, legend_y + 0.1, 'Match with prev', fontsize=8, va='center')

    ax.text(8.5, legend_y + 0.1, '0', fontsize=10, color=color_0, fontweight='bold', va='center')
    ax.text(9, legend_y + 0.1, '/', fontsize=10, va='center')
    ax.text(9.3, legend_y + 0.1, '1', fontsize=10, color=color_1, fontweight='bold', va='center')
    ax.text(9.8, legend_y + 0.1, 'alleles', fontsize=8, va='center')


def create_animation(haplotypes, output_file='pbwt_animation.gif', fps=1):
    """
    Create an animated GIF showing the PBWT algorithm frame by frame.
    """
    M, N = haplotypes.shape

    # Compute PBWT
    a_list, d_list = compute_pbwt(haplotypes)

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    fig.patch.set_facecolor('white')

    def animate(frame):
        k = frame
        create_figure1_frame(haplotypes, a_list, d_list, k, ax)
        return []

    # Create animation
    anim = FuncAnimation(
        fig, animate,
        frames=N,
        interval=1000 // fps,
        blit=False
    )

    # Save as GIF
    print(f"Creating animation with {N} frames...")
    writer = PillowWriter(fps=fps)
    anim.save(output_file, writer=writer, dpi=100)
    print(f"Animation saved to {output_file}")

    plt.close()
    return a_list, d_list


def save_individual_frames(haplotypes, output_dir='.', prefix='pbwt_frame'):
    """
    Save individual frames as PNG files.
    """
    import os

    M, N = haplotypes.shape

    # Compute PBWT
    a_list, d_list = compute_pbwt(haplotypes)

    print(f"Saving {N} individual frames...")

    for k in range(N):
        fig, ax = plt.subplots(figsize=(14, 8))
        fig.patch.set_facecolor('white')

        create_figure1_frame(haplotypes, a_list, d_list, k, ax)

        filename = os.path.join(output_dir, f'{prefix}_{k:02d}.png')
        fig.savefig(filename, dpi=150, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"  Saved {filename}")

    return a_list, d_list


def print_haplotype_panel(haplotypes):
    """Print the haplotype panel in a nice format."""
    M, N = haplotypes.shape
    print("\nHaplotype Panel:")
    print("-" * (N * 2 + 5))
    for i in range(M):
        hap_str = ' '.join(str(x) for x in haplotypes[i])
        print(f"h{i}: {hap_str}")
    print("-" * (N * 2 + 5))


def print_pbwt_state(haplotypes, a_list, d_list):
    """Print the PBWT state at each position."""
    M, N = haplotypes.shape

    print("\nPBWT Sort Order at Each Position:")
    print("=" * 50)

    for k in range(N + 1):
        a_k = a_list[k]
        d_k = d_list[k]

        print(f"\nPosition k={k}:")
        print(f"  Sort order a[k]: {list(a_k)}")
        print(f"  Divergence d[k]: {list(d_k)}")

        if k < N:
            print(f"  Sorted haplotypes at site {k}:")
            for i, hap_idx in enumerate(a_k):
                prefix = ''.join(str(x) for x in haplotypes[hap_idx, :k+1])
                print(f"    Row {i}: h{hap_idx} -> prefix: {prefix if prefix else '(empty)'}, "
                      f"allele at k={k}: {haplotypes[hap_idx, k]}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='PBWT Animation Generator')
    parser.add_argument('--haplotypes', '-n', type=int, default=8,
                        help='Number of haplotypes (default: 8)')
    parser.add_argument('--sites', '-s', type=int, default=12,
                        help='Number of sites (default: 12)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    parser.add_argument('--output', '-o', type=str, default='pbwt_animation.gif',
                        help='Output filename (default: pbwt_animation.gif)')
    parser.add_argument('--fps', type=int, default=1,
                        help='Frames per second (default: 1)')
    parser.add_argument('--frames-only', action='store_true',
                        help='Save individual frames instead of animation')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print detailed PBWT state information')

    args = parser.parse_args()

    # Generate haplotypes
    print(f"Generating {args.haplotypes} haplotypes with {args.sites} sites...")
    haplotypes = simulate_haplotypes(args.haplotypes, args.sites, args.seed)

    print_haplotype_panel(haplotypes)

    if args.frames_only:
        a_list, d_list = save_individual_frames(haplotypes, prefix='pbwt_frame')
    else:
        a_list, d_list = create_animation(haplotypes, args.output, args.fps)

    if args.verbose:
        print_pbwt_state(haplotypes, a_list, d_list)

    print("\nDone!")
