#!/usr/bin/env python
# coding: utf-8

import os
import re

import pandas as pd
from argparse import ArgumentParser, ArgumentTypeError


UNKNOWN_DIATOM = "Unknown diatom"
UNKNOWN_DIATOM_SOURCES = {
    "Bacillariophyceae_XXX",
    "Coscinodiscophyceae_XXX",
    "Mediophyceae_XXX",
}
# Same collapse as docs/replacement_taxa.tsv. Unknown diatom is included:
# those contigs are diatoms and stay annotated.
DIATOM_CLASS = "Diatoms"
DIATOM_CLASS_SOURCES = {
    "Bacillariophyceae",
    "Coscinodiscophyceae",
    "Mediophyceae",
    "Diatomeae_X",
    UNKNOWN_DIATOM,
    *UNKNOWN_DIATOM_SOURCES,
}
PLACEHOLDER_GENUS_RE = re.compile(r"_X+$")
PLASTID_GENUS_RE = re.compile(r":plas\b|:plastid\b", re.IGNORECASE)


def parse_bool(value):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n"}:
        return False
    raise ArgumentTypeError(f"expected true or false, got {value!r}")


def get_args():
    parser = ArgumentParser(description="Choose a subset of cells that pass a criteria")
    parser.add_argument("--path", help="Path for all the files", required=True)
    parser.add_argument(
        "--metapr2_file",
        help="metapr2_file = The output files for metapr2 blast file",
        required=True,
    )
    parser.add_argument(
        "--pr2_file",
        help="pr2_file =  The output files for pr2 blast file",
        required=True,
    )
    parser.add_argument(
        "--replacement_file",
        help="A tab delimited file for replacing taxonomic names and assigning colors for each taxonomic group",
        required=True,
    )
    parser.add_argument(
        "--identity",
        help="Threshold for read percent homology",
        default=99,
        type=float,
    )
    parser.add_argument(
        "--length",
        help="Threshold for alignment length",
        default=300,
        type=float,
    )
    parser.add_argument(
        "--level",
        help="Taxonomic level for the summary (e.g. Genus, Annotation)",
        default="Annotation",
    )
    parser.add_argument(
        "--force-annotation",
        type=parse_bool,
        default=True,
        help=(
            "If true (default), a real genus wins over a placeholder "
            "(*_X / *_XX / *_XXX / Unknown diatom) even at a lower bitscore, "
            "and placeholder contigs are dropped when the cell has a real genus. "
            "If false, the highest bitscore wins and only Unknown diatom "
            "co-hits are dropped."
        ),
    )
    return parser.parse_args()


def is_placeholder_taxon(name: object) -> bool:
    """True for unresolved PR2-style placeholders like Bacillariophyceae_XXX."""
    text = str(name)
    return bool(PLACEHOLDER_GENUS_RE.search(text)) or text == UNKNOWN_DIATOM


def load_data(path, file, db):
    list_columns = (
        "qseqid sseqid pident length mismatch gapopen "
        "qstart qend sstart send evalue bitscore"
    ).split(" ")

    print("reading file " + file + "...")
    full = file if os.path.isabs(file) else os.path.abspath(file)
    if os.path.getsize(full) == 0:
        print("empty blast file " + file)
        return pd.DataFrame(
            columns=list_columns
            + [
                "cell",
                "rrna_type",
                "rrna_loc",
                "Domain",
                "Clade",
                "Supergroup",
                "Division",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
                "db",
            ]
        )
    data = pd.read_table(full, sep="\t", names=list_columns)
    data["cell"] = data["qseqid"].str.split(r"i\d+_", expand=True)[1]
    data["cell"] = data["cell"].astype(str).str.replace(
        r":[0-9]+-[0-9]+\([+-]\)$", "", regex=True
    )
    if db == "PR2":
        data[
            [
                "sseqid",
                "rrna_type",
                "rrna_loc",
                "clone",
                "Domain",
                "Clade",
                "Supergroup",
                "Division",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
            ]
        ] = data["sseqid"].str.split("|", expand=True)
        # 18S only; exclude plastid 16S (and any non-18S). Nucleomorph 18S kept.
        before = len(data)
        data = data[
            (data["rrna_type"] == "18S_rRNA") & (data["rrna_loc"] != "plastid")
        ].copy()
        print(f"PR2 kept 18S non-plastid: {len(data)}/{before} rows")
    elif db == "metaPR2":
        data[
            [
                "sseqid",
                "Clade",
                "Supergroup",
                "Division",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species",
            ]
        ] = data["sseqid"].str.split("|", expand=True)
    else:
        print("Unrecognized data type! Exiting...")
        raise SystemExit(1)

    # Extra safety: drop plastid-tagged genus labels (e.g. Chaetoceros:plas).
    before = len(data)
    data = data[~data["Genus"].astype(str).str.contains(PLASTID_GENUS_RE, na=False)].copy()
    if before != len(data):
        print(f"Dropped plastid-tagged genera: {before - len(data)} rows")

    data["db"] = db
    print("Finished processing file")
    return data


def subset_results(sample_data):
    """Keep hits that pass length and identity. Genus choice is across both databases."""
    print("Filtering queries by length and identity...")
    sample_data = sample_data[
        (sample_data["length"] >= length) & (sample_data["pident"] >= identity)
    ]
    return sample_data


def prepare_db_hits(path, file, data_type, replacement_df, cells):
    """Load one BLAST table, apply replacements, and keep hits that pass the cutoffs."""
    data = load_data(path, file, data_type)
    global super_dict
    super_dict = replacement_process(data, replacement_df)
    if data.empty:
        return data
    data["Annotation"] = data["Family"].map(super_dict)
    data["barcode"] = data["cell"].str.split(".", expand=True)[0]
    if cells != -1:
        data = data[data["cell"].isin(cells)]
    data = subset_results(data)
    return data


def rename_unknown_diatom(df: pd.DataFrame) -> pd.DataFrame:
    """Map unresolved diatom class placeholders -> Unknown diatom."""
    out = df.copy()
    taxon = out[level].astype(str)
    out[level] = taxon.where(~taxon.isin(UNKNOWN_DIATOM_SOURCES), UNKNOWN_DIATOM)
    return out


def collapse_diatom_classes(df: pd.DataFrame) -> pd.DataFrame:
    """At Class, merge diatom classes, including Unknown diatom, into Diatoms."""
    if df.empty or level != "Class":
        return df
    out = df.copy()
    taxon = out[level].astype(str)
    n = int(taxon.isin(DIATOM_CLASS_SOURCES).sum())
    out[level] = taxon.where(~taxon.isin(DIATOM_CLASS_SOURCES), DIATOM_CLASS)
    print(f"Collapsed diatom classes to {DIATOM_CLASS}: {n} hit rows")
    return out


def drop_inconclusive_contigs(df: pd.DataFrame) -> pd.DataFrame:
    """June-style: drop contigs whose best bitscore hits disagree on taxonomy.

    If multiple distinct taxa share the maximum bitscore on a contig, the contig
    is inconclusive and is removed (not randomly assigned).
    """
    if df.empty:
        return df
    out = df.copy()
    out["_max_bits"] = out.groupby("qseqid")["bitscore"].transform("max")
    top = out[out["bitscore"] == out["_max_bits"]]
    n_taxa = top.groupby("qseqid")[level].nunique()
    inconclusive = set(n_taxa[n_taxa > 1].index)
    print(
        f"Dropped inconclusive contigs (multi-{level} at best bitscore): "
        f"{len(inconclusive)}"
    )
    kept = out[~out["qseqid"].isin(inconclusive)].drop(columns=["_max_bits"])
    return kept


def best_hit_per_contig(df: pd.DataFrame, force_annotation: bool) -> pd.DataFrame:
    """Keep one hit per contig across PR2∪metaPR2.

    With force annotation, a real genus that passed the filters wins even when
    a placeholder has a higher bitscore. Without it, the highest bitscore wins
    and a real genus is only a tie-break. Ties within the chosen group are
    broken by alignment length.
    """
    if df.empty:
        return df
    out = df.copy()
    out["_placeholder"] = out[level].map(is_placeholder_taxon).astype(int)
    if force_annotation:
        sort_by = ["qseqid", "_placeholder", "bitscore", "length"]
        ascending = [True, True, False, False]
    else:
        sort_by = ["qseqid", "bitscore", "length", "_placeholder"]
        ascending = [True, False, False, True]
    out = out.sort_values(by=sort_by, ascending=ascending)
    out = out.drop_duplicates(subset=["qseqid"], keep="first")
    return out.drop(columns=["_placeholder"])


def drop_unknown_diatom_when_other_taxa(df: pd.DataFrame) -> pd.DataFrame:
    """Drop Unknown diatom rows on cells that have any other taxon."""
    if df.empty:
        return df
    out = df.copy()
    is_unk = out[level].astype(str) == UNKNOWN_DIATOM
    has_named = out.assign(_named=~is_unk).groupby("cell")["_named"].transform("any")
    kept = out[~(is_unk & has_named)].copy()
    dropped = int(len(out) - len(kept))
    print(
        f"Dropped Unknown diatom co-hits on cells with another taxon: {dropped} contig rows"
    )
    return kept


def drop_placeholders_when_real_genus(df: pd.DataFrame) -> pd.DataFrame:
    """Drop placeholder contigs on cells that have a real genus.

    Cells whose only calls are placeholders are kept.
    """
    if df.empty:
        return df
    out = df.copy()
    is_ph = out[level].map(is_placeholder_taxon)
    has_real = out.assign(_real=~is_ph).groupby("cell")["_real"].transform("any")
    kept = out[~(is_ph & has_real)].copy()
    dropped = int(len(out) - len(kept))
    print(
        f"Dropped placeholder contigs on cells with a real genus: {dropped} contig rows"
    )
    return kept


def annotation_data(pr2_data, metapr2_data, force_annotation: bool):
    """Merge PR2 and metaPR2, then keep one hit per contig."""
    mode = "force annotation" if force_annotation else "best bitscore"
    print(f"Annotating cells (merge PR2+metaPR2, {mode}, one hit per contig)")
    cols = ["qseqid", "cell", level, "length", "bitscore", "pident", "db"]
    frames = []
    for data, name in ((pr2_data, "PR2"), (metapr2_data, "metaPR2")):
        if data is None or data.empty:
            print(f"{name}: empty after filters")
            continue
        sub = data.copy()
        missing = [c for c in cols if c not in sub.columns]
        if missing:
            raise SystemExit(f"ERROR: {name} missing columns {missing}")
        sub = sub[cols]
        sub["cell"] = sub["cell"].astype(str).str.strip()
        frames.append(sub)

    if not frames:
        return pd.DataFrame(columns=["cell", level, "percent", "length"])

    # 1) Merge like June: concatenate both databases.
    merged = pd.concat(frames, ignore_index=True)
    merged = rename_unknown_diatom(merged)
    # Class summary: Bacillariophyceae, Coscinodiscophyceae, Mediophyceae,
    # Diatomeae_X, and Unknown diatom are one annotated class.
    merged = collapse_diatom_classes(merged)

    if not force_annotation:
        # Best-bitscore mode: a contig whose top bitscore is shared by more
        # than one taxon is inconclusive and is dropped.
        merged = drop_inconclusive_contigs(merged)

    # One hit per contig. Force annotation lets a real genus beat a placeholder
    # even at a lower bitscore.
    merged = best_hit_per_contig(merged, force_annotation)

    if force_annotation:
        merged = drop_placeholders_when_real_genus(merged)
    else:
        merged = drop_unknown_diatom_when_other_taxa(merged)

    data_annotations = merged[["qseqid", "cell", level, "length"]].copy()
    data_annotations["percent"] = data_annotations["length"] / data_annotations.groupby(
        "cell"
    )["length"].transform("sum")
    data_annotations = (
        data_annotations.groupby(["cell", level], as_index=False)
        .agg({"percent": "sum", "length": "sum"})
    )
    return data_annotations


def replace(data, tax_level, name, new_name):
    group = set(data[data[tax_level] == name]["Family"])
    return {x: new_name for x in group}


def replacement_process(df_top, replacement_df):
    print("Start replacing process...")
    if df_top.empty or "Family" not in df_top.columns:
        print("Finished replacing process")
        return {}
    list_replacement = []
    for _index, row in replacement_df.iterrows():
        list_replacement.append(
            replace(df_top, row["level"], row["old_name"], row["new_name"])
        )

    super_dict = {key: val for d in list_replacement for key, val in d.items()}
    large_list = set(super_dict.keys())
    rest_species = {
        x: "Other eukaryotes"
        for x in [x for x in list(df_top.Family) if x not in large_list]
    }
    super_dict = {key: val for d in [super_dict, rest_species] for key, val in d.items()}
    print("Finished replacing process")
    return super_dict


def main():
    args = get_args()
    path = args.path
    metapr2_file = args.metapr2_file
    pr2_file = args.pr2_file
    replacement_file = args.replacement_file
    global level, length, identity
    length = args.length
    identity = args.identity
    level = args.level
    force_annotation = args.force_annotation
    cells = -1
    print(
        "Running blast summary for identity",
        identity,
        "and alignment length",
        length,
        "force_annotation",
        force_annotation,
    )
    if force_annotation:
        rule = (
            f"a real {level} wins over a placeholder even at a lower bitscore; "
            "placeholder contigs dropped when the cell has a real genus"
        )
    else:
        rule = (
            "highest bitscore wins (real genus only breaks ties); "
            "Unknown diatom dropped when the cell has another taxon"
        )
    print(
        "Merge mode: PR2+metaPR2 concat; one hit per contig; "
        f"{rule}; "
        f"{sorted(UNKNOWN_DIATOM_SOURCES)} -> {UNKNOWN_DIATOM}; "
        "PR2 18S only (non-plastid)"
        + (
            f"; {sorted(DIATOM_CLASS_SOURCES)} -> {DIATOM_CLASS}"
            if level == "Class"
            else ""
        )
    )
    if os.path.isabs(replacement_file):
        replacement_path = replacement_file
    else:
        replacement_path = os.path.abspath(replacement_file)
    names_df = pd.read_table(replacement_path, sep="\t")

    global metapr2_data, pr2_data
    metapr2_data = prepare_db_hits(
        path=path,
        file=metapr2_file,
        data_type="metaPR2",
        replacement_df=names_df,
        cells=cells,
    )
    pr2_data = prepare_db_hits(
        path=path,
        file=pr2_file,
        data_type="PR2",
        replacement_df=names_df,
        cells=cells,
    )

    data_annotations = annotation_data(pr2_data, metapr2_data, force_annotation)
    out_csv = os.path.join(path, "blast_summary." + level + ".csv")
    data_annotations.to_csv(out_csv, index=False)
    print(
        f"Wrote {out_csv} ({len(data_annotations)} rows, "
        f"{data_annotations['cell'].nunique()} cells)"
    )


if __name__ == "__main__":
    main()
