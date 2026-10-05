#!/usr/bin/env bash
# Set up a new dada2 analysis in the current directory.
# Usage (run in the project root):
#   curl -fsSL https://raw.githubusercontent.com/msuokas-bco/dada2_template/main/new_project.sh | bash
# Overwrite an existing template:
#   curl -fsSL https://raw.githubusercontent.com/msuokas-bco/dada2_template/main/new_project.sh | FORCE=1 bash

set -euo pipefail

url="https://raw.githubusercontent.com/msuokas-bco/dada2_template/main/dada2_template.qmd"
qmd="dada2_template.qmd"

if [[ -e "$qmd" && "${FORCE:-0}" != "1" ]]; then
    echo "[ERROR] $qmd already exists in $(pwd). Set FORCE=1 to overwrite." >&2
    exit 1
fi

# Download to a temp file first so a failed download never leaves a broken template
tmp=$(mktemp)
trap 'rm -f "$tmp"' EXIT
curl -fsSL "$url" -o "$tmp"
mv "$tmp" "$qmd"

mkdir -p reads results

version=$(grep -m1 '^version:' "$qmd" | sed 's/^version:[[:space:]]*//; s/"//g')

echo "Fetched $qmd (template version ${version:-unknown}) into $(pwd)"
echo
echo "Next steps:"
echo "  1. Copy metadata.tsv into the project root (needs a 'sampleid' column)"
echo "  2. Copy the demultiplexed FASTQ files into reads/"
echo "  3. quarto render $qmd -P amplicon:<target> [-P group:<metadata column>]"
echo "     (without group, the example analysis is skipped)"
