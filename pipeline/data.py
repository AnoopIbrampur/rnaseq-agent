"""Load the airway (GSE52778) raw counts, sample metadata, and gene-symbol map.

Biology notes
-------------
* A *count matrix* has one row per gene and one column per sample. Each cell
  is the number of sequencing reads that mapped to that gene in that sample,
  a rough measure of how much the gene was expressed.
* We keep the counts *raw* (integers). DESeq2 models counts directly with a
  negative binomial distribution and does its own normalization for
  sequencing depth internally. Feeding it pre-normalized values (FPKM, TPM)
  would break its statistical model, because the variance of a count depends
  on its absolute size, information that normalization throws away.
* NCBI's counts use NCBI Gene IDs (integers like 2289 for FKBP5). Pathway
  databases use gene *symbols* (FKBP5), so we map IDs to symbols with NCBI's
  own gene_info table.
"""

from __future__ import annotations

import gzip
import hashlib
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from pipeline import config


class DataIntegrityError(RuntimeError):
    """A raw input file is missing or does not match its recorded checksum."""


def sha256(path: Path) -> str:
    """Return the hex SHA-256 of a file, streamed in chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_raw_file(path: Path) -> str:
    """Check a raw file exists and matches config.EXPECTED_SHA256.

    Returns the checksum. Raises DataIntegrityError with a next step otherwise.
    """
    if not path.exists():
        raise DataIntegrityError(
            f"Missing raw file {path}. Download it as described in data/SOURCES.md "
            f"and place it in {config.DATA_RAW}/."
        )
    actual = sha256(path)
    expected = config.EXPECTED_SHA256.get(path.name)
    if expected is not None and actual != expected:
        raise DataIntegrityError(
            f"Checksum mismatch for {path.name}: expected {expected}, got {actual}. "
            "The file changed since it was recorded. Re-download it, or if the "
            "change is intended, update config.EXPECTED_SHA256 and data/SOURCES.md."
        )
    return actual


def parse_series_matrix(path: Path = config.SERIES_MATRIX_FILE) -> pd.DataFrame:
    """Parse sample metadata (all 16 samples) from a GEO series matrix file."""
    rows: dict[str, list[str]] = {}
    characteristics: list[list[str]] = []
    with gzip.open(path, "rt") as f:
        for line in f:
            if not line.startswith("!Sample_"):
                continue
            key, *values = line.rstrip("\n").split("\t")
            values = [v.strip('"') for v in values]
            if key == "!Sample_characteristics_ch1":
                characteristics.append(values)
            else:
                rows.setdefault(key, values)

    meta = pd.DataFrame(
        {"sample": rows["!Sample_geo_accession"], "title": rows["!Sample_title"]}
    )
    # Characteristics lines look like "cell line: N61311".
    for values in characteristics:
        name = values[0].split(":", 1)[0].strip().replace(" ", "_")
        meta[name] = [v.split(":", 1)[1].strip() for v in values]
    return meta.set_index("sample")


def load_metadata(path: Path = config.SERIES_MATRIX_FILE) -> pd.DataFrame:
    """Return metadata for the 8 dex/untreated samples.

    Columns: title, cell (donor line), dex ('untrt' or 'trt'), ercc_mix.
    `dex` is an ordered categorical with 'untrt' as the reference level, so a
    positive log2 fold change later means "higher with dexamethasone".
    """
    verify_raw_file(path)
    meta = parse_series_matrix(path)
    meta = meta[meta["treatment"].isin(config.KEEP_TREATMENTS)].copy()
    meta["dex"] = pd.Categorical(
        meta["treatment"].map(config.KEEP_TREATMENTS), categories=["untrt", "trt"]
    )
    meta = meta.rename(columns={"cell_line": "cell"})
    meta["cell"] = meta["cell"].astype("category")
    return meta[["title", "cell", "dex", "ercc_mix"]]


def load_counts(
    path: Path = config.COUNTS_FILE, samples: list[str] | None = None
) -> pd.DataFrame:
    """Return the raw integer count matrix (genes x samples), indexed by NCBI GeneID.

    If `samples` is given, columns are restricted to (and ordered as) that list.
    """
    verify_raw_file(path)
    counts = pd.read_csv(path, sep="\t", index_col="GeneID")
    counts.index = counts.index.astype(int)
    if samples is not None:
        missing = set(samples) - set(counts.columns)
        if missing:
            raise DataIntegrityError(
                f"Samples {sorted(missing)} are in the metadata but not the count "
                "matrix. The two files may come from different GEO releases."
            )
        counts = counts[samples]
    if (counts < 0).any().any() or not all(
        pd.api.types.is_integer_dtype(t) for t in counts.dtypes
    ):
        raise DataIntegrityError(
            f"{path.name} must contain non-negative integer counts. Check that you "
            "downloaded the raw_counts file, not FPKM/TPM."
        )
    return counts


def build_gene_map(
    gene_info_path: Path = config.GENE_INFO_FILE,
    gene_history_path: Path = config.GENE_HISTORY_FILE,
    out: Path = config.GENE_MAP_FILE,
) -> pd.DataFrame:
    """Build and cache an NCBI GeneID -> symbol (+ Ensembl ID) table.

    Current GeneIDs come from gene_info. The count matrix was annotated against
    an older NCBI release, so some of its GeneIDs have since been retired. For
    those, gene_history gives either the replacement GeneID (we take the
    current symbol) or the last symbol the retired gene had. `symbol_source`
    records which path each row took, so nothing is silently relabeled.
    """
    verify_raw_file(gene_info_path)
    verify_raw_file(gene_history_path)
    info = pd.read_csv(
        gene_info_path,
        sep="\t",
        usecols=["GeneID", "Symbol", "dbXrefs", "type_of_gene", "description"],
    )
    ensembl = info["dbXrefs"].str.extract(r"Ensembl:(ENSG\d+)", expand=False)
    current = pd.DataFrame(
        {
            "gene_id": info["GeneID"].astype(int),
            "symbol": info["Symbol"],
            "ensembl_id": ensembl,
            "gene_type": info["type_of_gene"],
            "description": info["description"],
            "symbol_source": "gene_info",
        }
    )

    hist = pd.read_csv(gene_history_path, sep="\t", dtype=str)
    hist = hist[~hist["Discontinued_GeneID"].astype(int).isin(current["gene_id"])]
    lookup = current.set_index("gene_id")
    retired = []
    for row in hist.itertuples(index=False):
        old_id = int(row.Discontinued_GeneID)
        if row.GeneID != "-" and int(row.GeneID) in lookup.index:
            new = lookup.loc[int(row.GeneID)]
            retired.append(
                (old_id, new.symbol, new.ensembl_id, new.gene_type, new.description,
                 f"gene_history:replaced_by_{row.GeneID}")
            )
        else:
            retired.append(
                (old_id, row.Discontinued_Symbol, None, None, None,
                 "gene_history:discontinued")
            )
    gene_map = pd.concat(
        [current, pd.DataFrame(retired, columns=current.columns)], ignore_index=True
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    gene_map.to_csv(out, index=False)
    return gene_map


def load_gene_map(path: Path = config.GENE_MAP_FILE) -> pd.DataFrame:
    """Load the cached gene map, building it from gene_info if absent."""
    if not path.exists():
        return build_gene_map(out=path)
    return pd.read_csv(path, dtype={"gene_id": int})


@dataclass
class Dataset:
    """The analysis-ready inputs: counts (genes x samples), metadata, gene map."""

    counts: pd.DataFrame
    metadata: pd.DataFrame
    gene_map: pd.DataFrame

    @property
    def input_checksums(self) -> dict[str, str]:
        """SHA-256 of each raw file this dataset was built from."""
        return {
            p.name: sha256(p)
            for p in (
                config.COUNTS_FILE,
                config.SERIES_MATRIX_FILE,
                config.GENE_INFO_FILE,
                config.GENE_HISTORY_FILE,
            )
        }


def load_dataset() -> Dataset:
    """Load counts for the 8 dex/untreated samples, their metadata, and the gene map."""
    metadata = load_metadata()
    counts = load_counts(samples=list(metadata.index))
    return Dataset(counts=counts, metadata=metadata, gene_map=load_gene_map())


def summarize(ds: Dataset) -> str:
    """Human-readable summary of the loaded dataset."""
    gm = ds.gene_map.set_index("gene_id")
    mapped = ds.counts.index.isin(gm.index)
    symbols = gm.reindex(ds.counts.index)["symbol"]
    dup = symbols[mapped].duplicated(keep=False).sum()
    nonzero = (ds.counts.sum(axis=1) > 0).sum()
    types = gm.reindex(ds.counts.index[mapped])["gene_type"].value_counts().head(5)
    sources = (
        gm.reindex(ds.counts.index[mapped])["symbol_source"]
        .str.replace(r"replaced_by_\d+", "replaced", regex=True)
        .value_counts()
    )
    lines = [
        "Sample metadata (8 samples kept of 16 in GEO):",
        ds.metadata.to_string(),
        "",
        f"Count matrix: {ds.counts.shape[0]:,} genes x {ds.counts.shape[1]} samples",
        f"Genes with any reads in these 8 samples: {nonzero:,}",
        f"Library sizes (total reads): "
        + ", ".join(f"{s}={n / 1e6:.1f}M" for s, n in ds.counts.sum().items()),
        f"GeneIDs mapped to a symbol: {mapped.sum():,} / {len(mapped):,}"
        f" ({mapped.mean():.2%}); duplicated symbols among mapped: {dup}",
        "Symbol source: " + ", ".join(f"{s}={n:,}" for s, n in sources.items()),
        "Top gene types: " + ", ".join(f"{t}={n:,}" for t, n in types.items()),
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(summarize(load_dataset()))
