#!/usr/bin/env python3
"""
Generate small synthetic MiSeq i100 amplicon datasets for testing the dada2 template.

Real amplicon sequences are extracted from a reference database (in-silico PCR with
the template primers) and turned into paired-end reads with i100 binned quality
scores (Q2/Q9/Q23/Q38) and quality-dependent substitution errors.

Datasets:
    16SV4  gut-like bacterial community from SILVA, 2x250 reads, fixed-length amplicons
    ITS2   soil/rhizosphere fungal community from UNITE, 2x300 reads, variable-length
           amplicons; amplicons shorter than the read length read through into the
           opposite primer and adapter

Output (relative to --outdir, default testdata/<region>):
    metadata.tsv                         sampleid + Sample_type
    reads/<sample>_S<n>_L001_R{1,2}_001.fastq.gz
    truth.tsv                            true amplicon (primer-free) counts per sample

Sequence sources:
    16SV4: SILVA SSU Ref NR99 138.2 (Quast et al. 2013, Nucleic Acids Res.
           41:D590-D596; https://www.arb-silva.de), licensed under CC BY 4.0.
    ITS2:  UNITE general FASTA release, dynamic SHs (Abarenkov et al. 2024, Nucleic
           Acids Res. 52:D791-D797; https://unite.ut.ee), licensed under CC BY-SA 4.0.
           UNITE sequences end at the ITS2/LSU boundary, so a conserved 5' LSU
           segment up to the ITS4 site is appended (Ascomycota / other phyla variant).

Only the Python standard library is needed. Fixed seed -> reproducible output.

Usage:
    python3 simulate_reads.py [--region 16SV4|ITS2] [--ref-fasta FILE] [--ref-taxa FILE]
                              [--outdir DIR] [--seed 2026]
"""

import argparse
import gzip
import os
import random
import re
from collections import Counter

# Primers and read lengths per region (primers as in the template registry)
REGIONS = {
    "16SV4":   {"fwd": "GTGYCAGCMGCCGCGGTAA",  "rev": "GGACTACNVGGGTWTCTAAT", "readlen": 250, "insert": (245, 260)},
    "16SV3V4": {"fwd": "CCTACGGGNGGCWGCAG",    "rev": "GACTACHVGGGTATCTAATCC", "readlen": 300, "insert": (395, 440)},
    "16SV4V5": {"fwd": "GTGYCAGCMGCCGCGGTAA",  "rev": "CCGTCAATTCMTTTGAGTTT",  "readlen": 300, "insert": (350, 390)},
    "ITS2":    {"fwd": "GTGARTCATCGARTCTTTG",  "rev": "TCCTCCGCTTATTGATATGC",  "readlen": 300, "insert": (150, 480)},
}

IUPAC = {"A": "A", "C": "C", "G": "G", "T": "T", "R": "AG", "Y": "CT", "S": "CG", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT"}
COMP = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")

# Illumina TruSeq adapter read-through (only used if the amplicon is shorter than the read)
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

# Conserved 5' LSU segment between the ITS2/LSU boundary and the ITS4 site, derived
# from UNITE sequences that extend into LSU
LSU_TAIL = {
    "p__Ascomycota": "GACCTCAAATCAGGTAGGATTACCCGCTGAACTTAA",
    "other":         "GACCTCAAATCAGGTAGGACTACCCGCTGAACTTAA",
}

DATASETS = {
    "16SV4": {
        "ref_fasta": "~/reference/silva-1382-seqs.fasta",
        "ref_taxa": "~/reference/silva-1382-taxa.tsv",
        # Taxa: (label, regex on the taxonomy string). Two Bacteroides species give
        # closely related ASVs; chloroplast and mitochondria test the non-target filter.
        "taxa": [
            ("Bacteroides_fragilis",     r"g__Bacteroides;s__Bacteroides_fragilis$"),
            ("Bacteroides_uniformis",    r"g__Bacteroides;s__Bacteroides_uniformis$"),
            ("Bacteroides_vulgatus",     r"g__Bacteroides;s__Bacteroides_vulgatus$"),
            ("Faecalibacterium",         r"g__Faecalibacterium;s__Faecalibacterium_prausnitzii$"),
            ("Blautia",                  r"g__Blautia;s__Blautia_obeum$"),
            ("Roseburia",                r"g__Roseburia;s__Roseburia_intestinalis$"),
            ("Akkermansia",              r"g__Akkermansia;s__Akkermansia_muciniphila$"),
            ("Alistipes",                r"g__Alistipes;s__Alistipes_putredinis$"),
            ("Bifidobacterium",          r"g__Bifidobacterium;s__Bifidobacterium_longum$"),
            ("Escherichia",              r"g__Escherichia-Shigella;s__Escherichia-Shigella_coli$"),
            ("Enterococcus",             r"g__Enterococcus;s__Enterococcus_faecalis$"),
            ("Streptococcus",            r"g__Streptococcus;s__Streptococcus_salivarius$"),
            ("Chloroplast",              r"o__Chloroplast;"),
            ("Mitochondria",             r"f__Mitochondria;"),
        ],
        # Relative abundance profiles. Treatment mimics an antibiotic-like disturbance.
        "profiles": {
            "Control": {
                "Bacteroides_fragilis": 10, "Bacteroides_uniformis": 12, "Bacteroides_vulgatus": 16,
                "Faecalibacterium": 20, "Blautia": 10, "Roseburia": 8, "Akkermansia": 6,
                "Alistipes": 7, "Bifidobacterium": 5, "Escherichia": 1.5, "Enterococcus": 0.5,
                "Streptococcus": 1.5, "Chloroplast": 1.5, "Mitochondria": 1,
            },
            "Treatment": {
                "Bacteroides_fragilis": 14, "Bacteroides_uniformis": 3, "Bacteroides_vulgatus": 6,
                "Faecalibacterium": 3, "Blautia": 4, "Roseburia": 1, "Akkermansia": 2,
                "Alistipes": 2, "Bifidobacterium": 10, "Escherichia": 28, "Enterococcus": 15,
                "Streptococcus": 10, "Chloroplast": 1.5, "Mitochondria": 1,
            },
        },
        "samples": [("Ctrl1", "Control"), ("Ctrl2", "Control"), ("Ctrl3", "Control"),
                    ("Trt1", "Treatment"), ("Trt2", "Treatment")],
        # A recurrent chimera (fixed parents and breakpoint) added to samples where both
        # parents are common, so that removeBimeraDenovo has a real ASV to flag.
        "recurrent_chimera": ("Faecalibacterium", "Bacteroides_vulgatus", 130, 4),
        # Dirichlet concentration: lower values give more sample-to-sample variation
        "concentration": 30,
    },
    "ITS2": {
        "ref_fasta": "~/reference/unite_dynamic.fasta",
        "ref_taxa": "~/reference/unite-taxa.tsv",
        # ITS2 lengths differ between taxa; two Penicillium species give close ASVs and
        # Malassezia is a low-level skin-derived contaminant
        "taxa": [
            ("Mortierella_alpina",          r"s__Mortierella_alpina$"),
            ("Solicoccozyma_terricola",     r"s__Solicoccozyma_terricola$"),
            ("Saitozyma_podzolica",         r"s__Saitozyma_podzolica$"),
            ("Penicillium_spinulosum",      r"s__Penicillium_spinulosum$"),
            ("Penicillium_canescens",       r"s__Penicillium_canescens$"),
            ("Pseudogymnoascus_roseus",     r"s__Pseudogymnoascus_roseus$"),
            ("Tetracladium",                r"g__Tetracladium;"),
            ("Fusarium_oxysporum",          r"s__Fusarium_oxysporum$"),
            ("Exophiala_salmonis",          r"s__Exophiala_salmonis$"),
            ("Trichoderma_harzianum",       r"s__Trichoderma_harzianum$"),
            ("Cladosporium_cladosporioides", r"s__Cladosporium_cladosporioides$"),
            ("Alternaria_alternata",        r"s__Alternaria_alternata$"),
            ("Russula_cyanoxantha",         r"s__Russula_cyanoxantha$"),
            ("Malassezia_restricta",        r"s__Malassezia_restricta$"),
        ],
        # Bulk soil dominated by saprotrophs and soil yeasts; rhizosphere enriched in
        # plant-associated fungi and an ectomycorrhizal Russula
        "profiles": {
            "Bulk_soil": {
                "Mortierella_alpina": 18, "Solicoccozyma_terricola": 14, "Saitozyma_podzolica": 10,
                "Penicillium_spinulosum": 9, "Penicillium_canescens": 6, "Pseudogymnoascus_roseus": 8,
                "Tetracladium": 8, "Fusarium_oxysporum": 4, "Exophiala_salmonis": 3,
                "Trichoderma_harzianum": 4, "Cladosporium_cladosporioides": 3,
                "Alternaria_alternata": 2, "Russula_cyanoxantha": 6, "Malassezia_restricta": 1,
            },
            "Rhizosphere": {
                "Mortierella_alpina": 8, "Solicoccozyma_terricola": 4, "Saitozyma_podzolica": 3,
                "Penicillium_spinulosum": 4, "Penicillium_canescens": 5, "Pseudogymnoascus_roseus": 2,
                "Tetracladium": 3, "Fusarium_oxysporum": 18, "Exophiala_salmonis": 12,
                "Trichoderma_harzianum": 12, "Cladosporium_cladosporioides": 8,
                "Alternaria_alternata": 6, "Russula_cyanoxantha": 12, "Malassezia_restricta": 1,
            },
        },
        "samples": [("Soil1", "Bulk_soil"), ("Soil2", "Bulk_soil"), ("Soil3", "Bulk_soil"),
                    ("Rhizo1", "Rhizosphere"), ("Rhizo2", "Rhizosphere")],
        "recurrent_chimera": ("Fusarium_oxysporum", "Exophiala_salmonis", 120, 4),
        "concentration": 80,
        "lsu_tail": LSU_TAIL,
    },
}

# MiSeq i100 bins (current software, introduced with the 1000-cycle kits): 2 (N only), 9, 23, 38.
# Approximate substitution rate per bin.
QBINS = {38: 2.0e-4, 23: 5.0e-3, 9: 1.0e-1}
N_RATE = 2.0e-5        # per-base no-call probability (Q2 'N')
CHIMERA_RATE = 0.02    # fraction of reads that are random-breakpoint chimeras
OFFTARGET_RATE = 0.01  # fraction of read pairs without primers (discarded by cutadapt)


def revcomp(s):
    return s.translate(COMP)[::-1]


def primer_regex(p):
    return "".join(f"[{IUPAC[b]}]" if len(IUPAC[b]) > 1 else b for b in p)


def resolve_primer(p, rng):
    """One concrete oligo from a degenerate primer pool."""
    return "".join(rng.choice(IUPAC[b]) for b in p)


def select_candidates(taxa_file, taxa, max_per_taxon=40):
    """Reference IDs per taxon, plus the phylum of each ID (for the LSU tail)."""
    patterns = [(lab, re.compile(rx)) for lab, rx in taxa]
    cands = {lab: [] for lab, _ in taxa}
    phylum = {}
    with open(taxa_file) as fh:
        next(fh)
        for line in fh:
            ref_id, tax = line.rstrip("\n").split("\t", 1)
            for lab, rx in patterns:
                if len(cands[lab]) < max_per_taxon and rx.search(tax):
                    cands[lab].append(ref_id)
                    phylum[ref_id] = tax.split(";")[1] if ";" in tax else ""
                    break
    return cands, phylum


def append_tail(seq, tail):
    """Append the LSU tail, merging any part of it already present at the 3' end."""
    for k in range(len(tail), 5, -1):
        if seq.endswith(tail[:k]):
            return seq + tail[k:]
    return seq + tail


def extract_amplicons(fasta, cands, phylum, fwd, rev, insert_range, lsu_tail=None):
    """In-silico PCR. For each taxon keep the most common error-free insert, or when
    every insert is unique, the one closest to the median length."""
    id2lab = {i: lab for lab, ids in cands.items() for i in ids}
    fwd_rx = re.compile(primer_regex(fwd))
    rev_rx = re.compile(primer_regex(revcomp(rev)))
    found = {lab: Counter() for lab in cands}

    def handle(rid, seq):
        lab = id2lab.get(rid)
        if lab is None:
            return
        seq = seq.upper().replace("U", "T")
        m1 = fwd_rx.search(seq)
        if not m1:
            return
        m2 = rev_rx.search(seq, m1.end())
        if m2:
            ins = seq[m1.end():m2.start()]
        elif lsu_tail:
            tail = lsu_tail.get(phylum.get(rid), lsu_tail["other"])
            ins = append_tail(seq[m1.end():], tail)
        else:
            return
        if insert_range[0] <= len(ins) <= insert_range[1] and set(ins) <= set("ACGT"):
            found[lab][ins] += 1

    rid, chunks = None, []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                if rid is not None:
                    handle(rid, "".join(chunks))
                rid, chunks = line[1:].split()[0], []
            elif rid in id2lab:
                chunks.append(line.strip())
        if rid is not None:
            handle(rid, "".join(chunks))

    amplicons = {}
    for lab, cnt in found.items():
        if not cnt:
            raise SystemExit(f"No amplicon found for {lab}; adjust the taxa patterns")
        best, n = cnt.most_common(1)[0]
        if n == 1 and len(cnt) > 1:
            lengths = sorted(len(s) for s in cnt)
            median = lengths[len(lengths) // 2]
            best = min(cnt, key=lambda s: (abs(len(s) - median), s))
        amplicons[lab] = best
    return amplicons


def quality_string(length, rng, read2):
    """Binned Q-scores: high Q38 throughout, Q23/Q9 more common toward the 3' end
    and in R2; low-quality calls tend to cluster; ~8% of reads are weaker overall."""
    factor = (1.45 if read2 else 1.0) * (2.5 if rng.random() < 0.08 else 1.0)
    quals, prev_low = [], False
    for i in range(length):
        f = i / length
        p9 = (0.002 + 0.010 * f ** 3) * factor
        p23 = (0.012 + 0.024 * f ** 2) * factor
        if prev_low:
            p9, p23 = p9 * 4, p23 * 3
        r = rng.random()
        q = 9 if r < p9 else 23 if r < p9 + p23 else 38
        quals.append(q)
        prev_low = q < 38
    return quals


def sequence_read(template, rng, read2):
    """Apply quality-dependent substitutions and rare no-calls."""
    quals = quality_string(len(template), rng, read2)
    bases = list(template)
    for i, q in enumerate(quals):
        if rng.random() < N_RATE:
            bases[i], quals[i] = "N", 2
        elif rng.random() < QBINS[q]:
            # transitions somewhat more common than transversions
            b = bases[i]
            trans = {"A": "G", "G": "A", "C": "T", "T": "C"}[b]
            others = [x for x in "ACGT" if x not in (b, trans)]
            bases[i] = trans if rng.random() < 0.5 else rng.choice(others)
    return "".join(bases), "".join(chr(q + 33) for q in quals)


def make_pair(insert, cfg, rng):
    fwd = resolve_primer(cfg["fwd"], rng)
    rev = resolve_primer(cfg["rev"], rng)
    L = cfg["readlen"]
    r1 = (fwd + insert + revcomp(rev) + ADAPTER_R1 + "A" * L)[:L]
    r2 = (rev + revcomp(insert) + revcomp(fwd) + ADAPTER_R2 + "A" * L)[:L]
    return r1, r2


def dirichlet(alpha, rng):
    g = [rng.gammavariate(a, 1) for a in alpha]
    s = sum(g)
    return [x / s for x in g]


def gzip_writer(path):
    # mtime=0 keeps regenerated files byte-identical
    return gzip.GzipFile(filename="", mode="wb", fileobj=open(path, "wb"), mtime=0)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--region", default="16SV4", choices=DATASETS)
    ap.add_argument("--ref-fasta", help="reference fasta (default depends on region)")
    ap.add_argument("--ref-taxa", help="reference taxonomy tsv (default depends on region)")
    ap.add_argument("--outdir", help="output directory (default testdata/<region>)")
    ap.add_argument("--seed", type=int, default=2026)
    args = ap.parse_args()

    ds = DATASETS[args.region]
    cfg = REGIONS[args.region]
    ref_fasta = os.path.expanduser(args.ref_fasta or ds["ref_fasta"])
    ref_taxa = os.path.expanduser(args.ref_taxa or ds["ref_taxa"])
    outdir = args.outdir or os.path.join(os.path.dirname(os.path.abspath(__file__)), args.region)

    # Separate streams: community composition stays fixed when the error model changes
    rng = random.Random(args.seed)
    srng = random.Random(args.seed + 1)

    print("Selecting reference sequences ...")
    cands, phylum = select_candidates(ref_taxa, ds["taxa"])
    amplicons = extract_amplicons(ref_fasta, cands, phylum, cfg["fwd"], cfg["rev"],
                                  cfg["insert"], ds.get("lsu_tail"))
    labels = [lab for lab, _ in ds["taxa"]]
    primers_len = len(cfg["fwd"]) + len(cfg["rev"])
    for lab in labels:
        amp = len(amplicons[lab]) + primers_len
        note = "  (read-through)" if amp < cfg["readlen"] else ""
        print(f"  {lab:30s} insert {len(amplicons[lab])} bp, amplicon {amp} bp{note}")

    readsdir = os.path.join(outdir, "reads")
    os.makedirs(readsdir, exist_ok=True)
    truth = {lab: {} for lab in labels}
    qstats = Counter()
    n_readthrough = 0

    for n, (sample, group) in enumerate(ds["samples"], start=1):
        depth = rng.randint(185, 215)
        # sample-level variation around the group profile (Dirichlet)
        base = ds["profiles"][group]
        tot = sum(base.values())
        props = dirichlet([ds["concentration"] * base[lab] / tot for lab in labels], rng)

        n_off = sum(rng.random() < OFFTARGET_RATE for _ in range(depth))
        n_chim = sum(rng.random() < CHIMERA_RATE for _ in range(depth - n_off))
        n_real = depth - n_off - n_chim
        draws = Counter(rng.choices(labels, weights=props, k=n_real))

        templates = []
        for lab in labels:
            truth[lab][sample] = draws[lab]
            templates += [amplicons[lab]] * draws[lab]

        # chimeras between two of the more abundant members
        top = [lab for lab, _ in draws.most_common(4)]
        for _ in range(n_chim):
            a, b = rng.sample(top, 2)
            bp = rng.randint(60, min(len(amplicons[a]), len(amplicons[b])) - 60)
            templates.append(amplicons[a][:bp] + amplicons[b][bp:])
        a, b, bp, copies = ds["recurrent_chimera"]
        if draws[a] >= 10 and draws[b] >= 10:
            templates += [amplicons[a][:bp] + amplicons[b][bp:]] * copies
            n_chim += copies
            depth += copies

        n_readthrough += sum(len(t) + primers_len < cfg["readlen"] for t in templates)
        pairs = [make_pair(t, cfg, srng) for t in templates]
        # off-target / primer-less pairs (random sequence)
        for _ in range(n_off):
            pairs.append(tuple("".join(srng.choice("ACGT") for _ in range(cfg["readlen"])) for _ in range(2)))
        srng.shuffle(pairs)

        f1 = os.path.join(readsdir, f"{sample}_S{n}_L001_R1_001.fastq.gz")
        f2 = os.path.join(readsdir, f"{sample}_S{n}_L001_R2_001.fastq.gz")
        with gzip_writer(f1) as o1, gzip_writer(f2) as o2:
            for k, (t1, t2) in enumerate(pairs, start=1):
                tile = srng.choice([1101, 1102, 1103, 1104, 2101, 2102, 2103, 2104])
                hdr = f"@SH00101:42:AAGTHK7M5:1:{tile}:{srng.randint(1000, 30000)}:{srng.randint(1000, 30000)}"
                s1, q1 = sequence_read(t1, srng, read2=False)
                s2, q2 = sequence_read(t2, srng, read2=True)
                qstats.update(q1 + q2)
                o1.write(f"{hdr} 1:N:0:{n}\n{s1}\n+\n{q1}\n".encode())
                o2.write(f"{hdr} 2:N:0:{n}\n{s2}\n+\n{q2}\n".encode())
        print(f"  {sample:6s} {group:11s} pairs={depth} (true={n_real}, chimeric={n_chim}, off-target={n_off})")

    with open(os.path.join(outdir, "metadata.tsv"), "w") as fh:
        fh.write("sampleid\tSample_type\n")
        for sample, group in ds["samples"]:
            fh.write(f"{sample}\t{group}\n")

    with open(os.path.join(outdir, "truth.tsv"), "w") as fh:
        fh.write("taxon\t" + "\t".join(s for s, _ in ds["samples"]) + "\tsequence\n")
        for lab in labels:
            fh.write(lab + "\t" + "\t".join(str(truth[lab][s]) for s, _ in ds["samples"])
                     + f"\t{amplicons[lab]}\n")

    total = sum(qstats.values())
    mean_q = sum((ord(c) - 33) * k for c, k in qstats.items()) / total
    q30 = sum(k for c, k in qstats.items() if ord(c) - 33 >= 30) / total
    print(f"Read pairs with primer/adapter read-through: {n_readthrough}")
    print(f"Mean Q {mean_q:.1f}, bases >= Q30 {100 * q30:.1f} %")
    print("Bin usage: " + ", ".join(f"Q{ord(c) - 33}: {100 * k / total:.2f} %"
                                    for c, k in sorted(qstats.items())))


if __name__ == "__main__":
    main()
