#!/usr/bin/env bash
# Set up the software environment for dada2_template.qmd.
# Usage:
#   curl -fsSL https://raw.githubusercontent.com/msuokas-bco/dada2_template/main/install.sh | bash
# Skip individual steps:
#   ... | SKIP_CONDA=1 SKIP_R=1 SKIP_REF=1 bash
#
# Steps:
#   1. conda environment 'python-311' in ~/.miniconda3 with cutadapt, mafft and FastTree
#      (Miniforge is installed into ~/.miniconda3 if it does not exist)
#   2. R packages from CRAN and Bioconductor (R itself must already be installed)
#   3. SILVA 138.2 reference files into ~/reference/ (UNITE must be downloaded manually)
# R and Quarto are not installed by this script; they are only checked.

set -euo pipefail

conda_root="$HOME/.miniconda3"
env_name="python-311"
ref_dir="$HOME/reference"

info() { echo "==> $*"; }
warn() { echo "[WARNING] $*" >&2; }
die()  { echo "[ERROR] $*" >&2; exit 1; }

command -v curl >/dev/null || die "curl is required."

# ---------------------------------------------------------------------------
# 1. conda environment with external tools
# ---------------------------------------------------------------------------
if [[ "${SKIP_CONDA:-0}" != "1" ]]; then
    if [[ ! -x "$conda_root/bin/conda" ]]; then
        [[ -e "$conda_root" ]] && die "$conda_root exists but contains no conda installation."
        info "Installing Miniforge into $conda_root"
        installer=$(mktemp)
        trap 'rm -f "$installer"' EXIT
        curl -fsSL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname -s)-$(uname -m).sh" \
            -o "$installer"
        bash "$installer" -b -p "$conda_root"
        rm -f "$installer"
    fi

    conda="$conda_root/bin/conda"
    # conda-forge and bioconda only; avoids the Anaconda defaults channel
    channels=(--override-channels -c conda-forge -c bioconda)
    tools=(cutadapt mafft fasttree)

    if [[ -d "$conda_root/envs/$env_name" ]]; then
        info "Updating conda environment '$env_name'"
        "$conda" install -y -n "$env_name" "${channels[@]}" "${tools[@]}"
    else
        info "Creating conda environment '$env_name'"
        "$conda" create -y -n "$env_name" "${channels[@]}" python=3.11 "${tools[@]}"
    fi

    env_bin="$conda_root/envs/$env_name/bin"
    "$env_bin/cutadapt" --version >/dev/null || die "cutadapt is not working in '$env_name'."
    [[ -x "$env_bin/mafft" ]] || die "mafft is missing from '$env_name'."
    [[ -x "$env_bin/FastTree" || -x "$env_bin/FastTreeMP" ]] || die "FastTree is missing from '$env_name'."
    info "cutadapt $("$env_bin/cutadapt" --version), mafft and FastTree are available in '$env_name'"
fi

# ---------------------------------------------------------------------------
# 2. R packages
# ---------------------------------------------------------------------------
if [[ "${SKIP_R:-0}" != "1" ]]; then
    command -v Rscript >/dev/null || die "R is not installed. Install it from https://cran.r-project.org and rerun (or set SKIP_R=1)."
    info "Installing R packages (already installed packages are skipped)"
    Rscript --vanilla - <<'EOF'
repos <- c(CRAN = "https://cloud.r-project.org")
options(repos = repos, Ncpus = max(1L, parallel::detectCores() - 1L))
lib <- Sys.getenv("R_LIBS_USER")
if (file.access(.libPaths()[1], 2) != 0) {
  dir.create(lib, recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(lib, .libPaths()))
}
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

pkgs <- c(
  # Bioconductor
  "dada2", "mia", "TreeSummarizedExperiment", "scater", "Biostrings",
  # CRAN
  "tidyverse", "vegan", "ape", "phangorn", "kableExtra", "patchwork",
  "ggrepel", "ggpubr", "ggsci", "bslib",
  "sessioninfo", "rmarkdown", "knitr"
)
missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  BiocManager::install(missing, ask = FALSE, update = FALSE)
}
failed <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(failed) > 0) {
  stop("These R packages failed to install: ", paste(failed, collapse = ", "),
       "\nOn Linux, check that system libraries (libcurl, libssl, libxml2, ",
       "libfontconfig, libharfbuzz, libfribidi, libfreetype, libpng, libtiff, libjpeg) are installed.")
}
message("All R packages are installed.")
EOF
fi

# ---------------------------------------------------------------------------
# 3. Reference databases
# ---------------------------------------------------------------------------
if [[ "${SKIP_REF:-0}" != "1" ]]; then
    mkdir -p "$ref_dir"
    # SILVA 138.2 formatted for DADA2 (https://doi.org/10.5281/zenodo.14169026)
    silva_url="https://zenodo.org/records/14169026/files"
    silva_files=(
        "silva_nr99_v138.2_toGenus_trainset.fa.gz 1764e2a36b4500ccb1c7d5261948a414"
        "silva_v138.2_assignSpecies.fa.gz 62fd939bd1aa9832d2089847f9f0b87d"
    )
    for entry in "${silva_files[@]}"; do
        read -r file md5 <<< "$entry"
        if [[ -s "$ref_dir/$file" ]]; then
            info "$file already exists in $ref_dir"
            continue
        fi
        info "Downloading $file"
        curl -fL --retry 3 "$silva_url/$file?download=1" -o "$ref_dir/$file.part"
        if command -v md5sum >/dev/null; then
            got=$(md5sum "$ref_dir/$file.part" | cut -d' ' -f1)
        else
            got=$(md5 -q "$ref_dir/$file.part")
        fi
        [[ "$got" == "$md5" ]] || { rm -f "$ref_dir/$file.part"; die "Checksum mismatch for $file."; }
        mv "$ref_dir/$file.part" "$ref_dir/$file"
    done

    unite="sh_general_release_dynamic_all_19.02.2025.fasta"
    if [[ ! -s "$ref_dir/$unite" ]]; then
        warn "UNITE reference for ITS targets is missing and cannot be downloaded automatically."
        warn "Download the UNITE general FASTA release (all eukaryotes, 19.02.2025) from"
        warn "https://unite.ut.ee/repository.php and place $unite in $ref_dir/"
    fi
fi

# ---------------------------------------------------------------------------
# Quarto check
# ---------------------------------------------------------------------------
if command -v quarto >/dev/null; then
    info "Quarto $(quarto --version) found"
else
    warn "Quarto is not installed. Install it from https://quarto.org/docs/get-started/"
fi

echo
info "Done. Start a new project with:"
echo "  curl -fsSL https://raw.githubusercontent.com/msuokas-bco/dada2_template/main/new_project.sh | bash"
