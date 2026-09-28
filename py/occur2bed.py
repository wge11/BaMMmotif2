"""Convert a BaMM!motif .occurrence file into a BED-like file.

The sequence headers of the scanned FASTA file must contain genomic
coordinates in the ``chrom:start-end`` form written by ``bedtools getfasta``
(0-based start, exclusive end), for example ``>chr1:1000-1205``.

Output columns (tab-separated, with a header line):

    #CHROM  START  END  STRAND  Pval

START is 0-based and END is exclusive, as in BED. Occurrences on the reverse
strand are mapped back onto forward-strand genomic coordinates.

Usage:
    python occur2bed.py SAMPLE_motif_1.occurrence [-o OUTDIR]
"""

import argparse
import os
import re
import sys

HEADER_RE = re.compile(r"^>?(?P<chrom>[^:\s]+):(?P<start>\d+)-(?P<end>\d+)")


def create_parser():
    parser = argparse.ArgumentParser(
        description="Convert a BaMM!motif .occurrence file into a BED-like file."
    )
    parser.add_argument(
        "occurrence_file", help=".occurrence file written by BaMMScan or BaMMmotif --scoreSeqset"
    )
    parser.add_argument(
        "-o", dest="outdir", default=None, help="output directory (default: next to the input file)"
    )
    return parser


def occurrence_to_bed(header, seq_length, strand, pos, pvalue):
    """Map one occurrence to (chrom, start, end, strand, pvalue).

    ``pos`` is the ``start..end`` field of the .occurrence file. Positions are
    1-based and inclusive on the sequence that BaMM!motif scans: the forward
    strand (positions 1..L), a separator, and the reverse complement
    (positions L+2..2L+1).
    """
    match = HEADER_RE.match(header)
    if match is None:
        raise ValueError(
            f"sequence header {header!r} contains no genomic coordinates (expected e.g. '>chr1:1000-1205')"
        )
    chrom = match.group("chrom")
    region_start = int(match.group("start"))
    motif_start, motif_end = (int(x) for x in pos.split(".."))

    if strand == "+":
        start = region_start + motif_start - 1
        end = region_start + motif_end
    else:
        # reverse-complement position p (1-based) corresponds to forward
        # position 2L + 2 - p (1-based)
        start = region_start + 2 * seq_length + 1 - motif_end
        end = region_start + 2 * seq_length + 2 - motif_start
    return chrom, start, end, strand, pvalue


def convert(occurrence_file, bed_file):
    with open(occurrence_file) as fin, open(bed_file, "w") as fout:
        print("#CHROM", "START", "END", "STRAND", "Pval", sep="\t", file=fout)
        header_line = fin.readline()
        if not header_line.startswith("seq\t"):
            raise ValueError(f"{occurrence_file} does not look like a .occurrence file")
        for line in fin:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            header, length, strand, pos, _pattern, pvalue = fields[:6]
            chrom, start, end, strand, pvalue = occurrence_to_bed(header, int(length), strand, pos, pvalue)
            print(chrom, start, end, strand, f"{float(pvalue):0.2e}", sep="\t", file=fout)


def main(argv=None):
    args = create_parser().parse_args(argv)

    outdir = args.outdir if args.outdir is not None else os.path.dirname(args.occurrence_file)
    if outdir and not os.path.exists(outdir):
        os.makedirs(outdir)
    basename = os.path.splitext(os.path.basename(args.occurrence_file))[0]
    bed_file = os.path.join(outdir, basename + ".bed")

    try:
        convert(args.occurrence_file, bed_file)
    except ValueError as err:
        sys.exit(f"Error: {err}")


if __name__ == "__main__":
    main()
