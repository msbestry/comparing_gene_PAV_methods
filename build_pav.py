# /// script
# requires-python = ">=3.9"
# dependencies = [
#     "pandas>=2.0",
# ]
# ///
"""

Run with uv: uv run scripts/mitch_pav_tables/build_pav_matrix.py

Inputs  (in data/mitch_pav_tables/):
    minigraph_pav.tsv, orthofinder_pav.tsv, pggb_pav.tsv

Final output:
    master_pav_matrix.tsv
"""

from pathlib import Path
import pandas as pd

# Inputs.
REPO_ROOT = Path(__file__).resolve().parent.parent.parent
DATA_DIR = REPO_ROOT / "data"
PAV_DIR = DATA_DIR / "mitch_pav_tables"

# Accession naming
CANONICAL_ACCESSIONS = ["a55015", "Bristol", "DSV-1", "DSV-2", "ERA3543"]

# Different sources spell the same accession in different ways so map every variant seen in the raw tables onto a
# standardized name.
ACCESSION_ALIASES = {
    "a55015": "a55015",
    "Bristol": "Bristol",
    "DSV-1": "DSV-1",
    "DSV_1": "DSV-1",
    "DSV_1_20210603": "DSV-1",
    "DSV-2": "DSV-2",
    "DSV_2": "DSV-2",
    "ERA3543": "ERA3543",
}

# Define Reference to anchor calls.
TARGET_ACCESSION = "a55015"

# Normalise funcs.
def canonical_accession(label):
    return ACCESSION_ALIASES.get(label, label)

def origin_from_seqid(seqid):
    """Extract the origin accession from a minigraph AnnotSeqID.
    The ID is '<accession>_<scaffold-number>'.
    """
    return canonical_accession(seqid.rsplit("_", 1)[0])

def origin_from_piped(name):
    """Extract the origin accession from a '<accession>|<gene-id>' name."""
    return canonical_accession(name.split("|", 1)[0])

def gene_token_from_piped(name):
    """Extract the gene-id portion from a '<accession>|<gene-id>' name."""
    return name.split("|", 1)[1] if "|" in name else name

def build_output(df, gene_ids, origins, pav_columns):
    """Assemble the common-schema table for one source.

    - Builds the ``gene`` key as '<origin>__<gene-id>'.
    - Copies the PAV flag columns into canonical accession order.
    - Keeps only rows whose origin is the TARGET_ACCESSION.
    """
    origins = pd.Series(list(origins), index=df.index)
    gene_ids = pd.Series(list(gene_ids), index=df.index)

    out = pd.DataFrame(index=df.index)
    out["gene"] = origins.str.cat(gene_ids, sep="__")

    raw_for_canon = {canon: raw for raw, canon in pav_columns.items()}
    for canon in CANONICAL_ACCESSIONS:
        out[canon] = df[raw_for_canon[canon]].values

    out = out[origins.values == TARGET_ACCESSION].reset_index(drop=True)
    return out

# These functions contain standardisation logic for the three inputs.
def process_minigraph(path):
    """Normalise the minigraph PAV table.

    Layout: first 4 columns are metadata (incl. GeneID and AnnotSeqID),
    remaining columns are per-accession PAV flags. Origin comes from AnnotSeqID.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    pav_columns = {c: canonical_accession(c) for c in df.columns[4:]}
    origins = df["AnnotSeqID"].map(origin_from_seqid)
    return build_output(df, df["GeneID"], origins, pav_columns)

def process_orthofinder(path):
    """Normalise the orthofinder PAV table.

    Layout: first column ``gene`` is '<accession>|<gene-id>', remaining columns
    are per-accession PAV flags.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    pav_columns = {c: canonical_accession(c) for c in df.columns[1:]}
    origins = df["gene"].map(origin_from_piped)
    gene_ids = df["gene"].map(gene_token_from_piped)
    return build_output(df, gene_ids, origins, pav_columns)

def process_pggb(path):
    """Normalise the pggb PAV table.

    Layout: first 4 columns are metadata (incl. ``name`` as '<accession>|<gene-id>'),
    remaining columns are per-accession PAV flags.
    """
    df = pd.read_csv(path, sep="\t", dtype=str)
    pav_columns = {c: canonical_accession(c) for c in df.columns[4:]}
    origins = df["name"].map(origin_from_piped)
    gene_ids = df["name"].map(gene_token_from_piped)
    return build_output(df, gene_ids, origins, pav_columns)


# Maps each raw input file to the function that knows how to normalize it.
PROCESSORS = {
    "minigraph_pav.tsv": process_minigraph,
    "orthofinder_pav.tsv": process_orthofinder,
    "pggb_pav.tsv": process_pggb,
}

# Final merge.
def normalise_tables():
    """Rewrite every raw PAV table into the common schema.

    Returns a dict mapping source prefix -> normalised output filename, for the
    merge step to consume. Missing input files are skipped with a warning.
    """
    print(f"Normalising PAV tables in {PAV_DIR}")
    normalised = {}
    for filename, processor in PROCESSORS.items():
        in_path = PAV_DIR / filename
        if not in_path.exists():
            print(f"  [skip] {filename} not found")
            continue

        out = processor(in_path)
        out_path = PAV_DIR / f"{in_path.stem}.normalised{in_path.suffix}"
        out.to_csv(out_path, sep="\t", index=False)
        print(f"  [ok]   {filename}: {len(out):>7,} {TARGET_ACCESSION} genes -> {out_path.name}")

        # Source prefix is the leading token of the filename (e.g. 'minigraph').
        prefix = filename.split("_", 1)[0]
        normalised[prefix] = out_path.name

    return normalised

# Merge to make one output.
def build_master_matrix(sources):
    """Outer-join the normalised source tables into one master matrix.

    ``sources`` maps source prefix -> normalised filename. Each source's PAV
    columns are prefixed with the source name to keep them distinct, and genes
    missing from a source are filled with 0 (absent).
    """
    print("\nBuilding master matrix")
    master = None
    for prefix, filename in sources.items():
        df = pd.read_csv(PAV_DIR / filename, sep="\t", dtype=str)
        # Drop any duplicate gene ids so the join stays one row per gene.
        df = df.drop_duplicates(subset="gene")
        pav_cols = [c for c in df.columns if c != "gene"]
        df = df.rename(columns={c: f"{prefix}_{c}" for c in pav_cols})

        master = df if master is None else master.merge(df, on="gene", how="outer")

    if master is None:
        print("  [warn] no normalised tables to merge; nothing written")
        return

    # Genes missing from a source merge in as NaN; treat absence as 0.
    master = master.fillna("0")

    out_path = PAV_DIR / "master_pav_matrix.tsv"
    master.to_csv(out_path, sep="\t", index=False)
    print(f"  [ok] {len(master):,} genes x {master.shape[1]} columns -> {out_path.name}")
    print(f"       columns: {list(master.columns)}")

    # Per-column summary of 1 / 0 counts across all genes (nulls are now 0).
    pav_cols = [c for c in master.columns if c != "gene"]
    print(f"\n  {'column':<40} {'ones':>8} {'zeros':>8}")
    for col in pav_cols:
        counts = master[col].value_counts()
        ones = int(counts.get("1", 0))
        zeros = int(counts.get("0", 0))
        print(f"  {col:<40} {ones:>8,} {zeros:>8,}")

# Define main()
def main():
    sources = normalise_tables()
    build_master_matrix(sources)

if __name__ == "__main__":
    main()