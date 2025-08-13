#!/usr/bin/env python3
"""
Utilities to parse bcftools AF stats and build a folded SFS (no plotting).
"""

from typing import Iterable, List, Tuple, Dict
from collections import defaultdict
import argparse
import sys


def parse_bcftools_af_data(data_text: str) -> List[Tuple[float, int]]:
    """
    Parse bcftools stats lines that start with 'AF' and extract allele
    frequency and count.

    Expected line format (tab-delimited):
    AF  <bin>  <allele_frequency>  <count>  ...

    Returns a list of tuples: (allele_frequency, count)
    """
    lines = data_text.strip().split("\n")
    af_items: List[Tuple[float, int]] = []
    for line in lines:
        if line.startswith("AF\t") and not line.startswith("# AF"):
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            try:
                af = float(parts[2])
                cnt = int(parts[3])
            except ValueError:
                continue
            af_items.append((af, cnt))
    return af_items


def fold_sfs(
    af_data: Iterable[Tuple[float, int]],
    num_bins: int = 20,
) -> Dict[int, int]:
    """
    Create a folded site frequency spectrum.

    - af_data: iterable of (allele_frequency, count)
    - num_bins: number of bins over 0..0.5 (inclusive of 0, exclusive of
      0.5 upper edge)

    Returns a dict mapping bin index [0..num_bins-1] -> count
    """
    folded_counts: Dict[int, int] = defaultdict(int)
    for af, cnt in af_data:
        # Fold frequencies > 0.5
        folded_af = 1.0 - af if af > 0.5 else af
        # Clamp to [0, 0.5)
        if folded_af < 0:
            folded_af = 0.0
        if folded_af >= 0.5:
            folded_af = 0.4999999999
        # Bin (*2 because original range is 0..1)
        bin_idx = int(folded_af * num_bins * 2)
        if bin_idx >= num_bins:
            bin_idx = num_bins - 1
        folded_counts[bin_idx] += cnt
    return folded_counts


def create_sfs_arrays(
    folded_counts: Dict[int, int],
    num_bins: int = 20,
) -> Tuple[List[float], List[int]]:
    """Convert folded counts to frequency centers and counts arrays."""
    frequencies: List[float] = []
    counts: List[int] = []
    # Bin width
    bin_width = 0.5 / num_bins
    for i in range(num_bins):
        freq = (i + 0.5) * bin_width
        cnt = folded_counts.get(i, 0)
        frequencies.append(freq)
        counts.append(cnt)
    return frequencies, counts


def process_full_bcftools_file(
    filename: str,
    num_bins: int = 20,
) -> Tuple[List[float], List[int]]:
    """Process a full bcftools stats file and return folded SFS arrays."""
    with open(filename, "r", encoding="utf-8") as f:
        content = f.read()
    af_items = parse_bcftools_af_data(content)
    folded_counts = fold_sfs(af_items, num_bins=num_bins)
    return create_sfs_arrays(folded_counts, num_bins=num_bins)


def save_folded_sfs(
    frequencies: List[float],
    counts: List[int],
    output_file: str,
) -> None:
    """Save folded SFS to a tab-delimited file."""
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("Folded_AF\tCount\n")
        for freq, cnt in zip(frequencies, counts):
            f.write(f"{freq:.6f}\t{int(cnt)}\n")


def main(argv=None) -> int:
    """CLI entry point: fold SFS from a bcftools stats file or demo input."""
    parser = argparse.ArgumentParser(
        description=(
            "Fold bcftools AF SFS from stats output (no plotting)"
        )
    )
    parser.add_argument(
        "input",
        nargs="?",
        help=(
            "Path to bcftools stats file. If omitted, runs a small demo."
        ),
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--bins",
        type=int,
        help="Number of bins for 0..0.5 (e.g., 20)",
    )
    group.add_argument(
        "--bin-width",
        dest="bin_width",
        type=float,
        help=(
            "Bin width over 0..0.5 (e.g., 0.025 => 20 bins)."
        ),
    )
    parser.add_argument(
        "--out",
        help="Optional path to save folded SFS as TSV",
    )
    args = parser.parse_args(argv)

    # Determine effective number of bins
    if args.bin_width is not None:
        if args.bin_width <= 0:
            print("--bin-width must be > 0", file=sys.stderr)
            return 2
        bins = max(1, int(0.5 / args.bin_width))
    else:
        bins = args.bins if args.bins is not None else 20

    if args.input:
        freqs, cnts = process_full_bcftools_file(
            args.input, num_bins=bins
        )
        print("Folded Site Frequency Spectrum (first 10 bins):")
        print("Frequency\tCount")
        for i in range(min(10, len(freqs))):
            print(f"{freqs[i]:.4f}\t\t{int(cnts[i])}")
        if args.out:
            save_folded_sfs(freqs, cnts, args.out)
    else:
        # Demo with a tiny inline sample
        sample = (
            "AF\t0\t0.000000\t5574\t3140\t2434\t0\t0\t0\t0\n"
            "AF\t0\t0.003559\t19097\t11060\t8037\t0\t0\t0\t0\n"
            "AF\t0\t0.007117\t11337\t6649\t4688\t0\t0\t0\t0\n"
            "AF\t0\t0.010676\t7354\t4443\t2911\t0\t0\t0\t0\n"
        )
        items = parse_bcftools_af_data(sample)
        folded = fold_sfs(items, num_bins=bins)
        freqs, cnts = create_sfs_arrays(folded, num_bins=bins)
        print("Parsed AF data (first 4 entries):")
        for af, c in items:
            print(f"AF: {af:.6f}, Count: {c}")
        print("\nFolded Site Frequency Spectrum (first 10 bins):")
        print("Frequency\tCount")
        for i in range(min(10, len(freqs))):
            print(f"{freqs[i]:.4f}\t\t{int(cnts[i])}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
