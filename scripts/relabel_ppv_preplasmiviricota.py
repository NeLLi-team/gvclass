#!/usr/bin/env python3
"""GVClass PPV rename: unify the Preplasmiviricota domain prefix PLV+VP -> PPV.

After issue17 (VP__->PLV__) + collision-resolution + ICTV-2025, the Preplasmiviricota
phylum spans two classes (Polintoviricetes, Virophaviricetes) under ONE domain prefix
that is still the class-level misnomer "PLV" (Polinton-like virus). This renames the
unified domain prefix to "PPV" (PreplasmiViricota) across all DB artifacts so the domain
label reflects the phylum, not one of its classes.

GVClass derives the domain from the reference ID prefix before "__" (summarize_full.py),
so IDs and tree leaves MUST be renamed, not just the lineage column.

Rename rule (applied per "|"-delimited token / per ID), double-underscore IDs and the
leading domain token ONLY:

    token == "PLV" or "VP"      -> "PPV"            (leading lineage domain token)
    token startswith "PLV__"    -> "PPV__" + rest   (IDs, tree leaves, pos7 organism)
    token startswith "VP__"     -> "PPV__" + rest

KEEP unchanged (NOT the domain):
  * marker HMM names: PLV_MCP_*, plv_mcp_caps_*, VP_MCP_*, VP_ATPase_*, VP_Penton_*,
    VP_PRO_*, gamadvirusMCP  (single-underscore; combined.hmm + faa filenames)
  * intra-lineage lifestyle tags: PLV_unclassified, VP_unclassified (single-underscore)
  * class names Polintoviricetes / Virophaviricetes and every other taxon.

The rule is exact: single-underscore tokens never start with "PLV__"/"VP__", so marker
names and lifestyle tags are provably untouched. "__vpdup" suffixes survive (prefix-only
replace).

Targets (resources/):
  labels.tsv               (runtime)    col0 id, col1 lineage
  inactive_labels.tsv      (archival)   col0 id, col1 lineage, col2 reason  [CRLF preserved]
  aliases.tsv              (build)      col0/1 ids, col3/4 lineages
  label_context.tsv        (build)      col0 id, col2/3 generic id-or-lineage
  faa_rewrite_map.tsv      (provenance) col1/2 headers, col4/5 ids
  database/faa/*.faa        (runtime)   ">" headers

Provenance note: historical-provenance domain tokens (faa_rewrite_map old_header,
label_context original_*_taxonomy values) are INTENTIONALLY migrated too, so the live DB
holds zero PLV/VP domain tokens anywhere (a clean, checkable global invariant). __vpdup and
all non-domain content are preserved.

Safety: writes are atomic (write *.tmp_ppv then os.replace); each modified file is backed up
to *.bak_ppv before its first write; original line terminators are preserved. Idempotent: if
PPV is already present this is a no-op. To guard the backup<->live invariant, --apply refuses
to run if any *.bak_ppv already exists (stale/partial prior run) unless --force is given.

Usage:
    python scripts/relabel_ppv_preplasmiviricota.py            # dry-run (no writes)
    python scripts/relabel_ppv_preplasmiviricota.py --apply    # write; back up to *.bak_ppv
    python scripts/relabel_ppv_preplasmiviricota.py --apply --force   # ignore stale *.bak_ppv
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "resources"
APPLY = "--apply" in sys.argv
FORCE = "--force" in sys.argv
BAK = ".bak_ppv"
TMP = ".tmp_ppv"


# --- token transforms --------------------------------------------------------
def remap_token(tok: str) -> str:
    if tok == "PLV" or tok == "VP":
        return "PPV"
    if tok.startswith("PLV__"):
        return "PPV__" + tok[5:]
    if tok.startswith("VP__"):
        return "PPV__" + tok[4:]
    return tok


def remap_lineage(lin: str) -> str:
    return "|".join(remap_token(f) for f in lin.split("|"))


def remap_header(hdr: str) -> str:
    # faa header: genome-id is the part before the first "|"; protein suffix untouched.
    if "|" in hdr:
        idpart, prot = hdr.split("|", 1)
        return remap_token(idpart) + "|" + prot
    return remap_token(hdr)


def remap_field(f: str) -> str:
    # generic: a lineage if it has "|", else a single token.
    return remap_lineage(f) if "|" in f else remap_token(f)


# --- per-line row transforms (return new line, no terminator) ----------------
def row_labels(line: str) -> str:
    p = line.split("\t")
    p[0] = remap_token(p[0])
    if len(p) > 1:
        p[1] = remap_lineage(p[1])
    return "\t".join(p)


row_inactive = row_labels  # same shape for the first two columns


def row_aliases(line: str) -> str:
    p = line.split("\t")
    if p:
        p[0] = remap_token(p[0])
    if len(p) > 1:
        p[1] = remap_token(p[1])
    if len(p) > 3:
        p[3] = remap_lineage(p[3])
    if len(p) > 4:
        p[4] = remap_lineage(p[4])
    return "\t".join(p)


def row_label_context(line: str) -> str:
    p = line.split("\t")
    if p:
        p[0] = remap_token(p[0])
    if len(p) > 2:
        p[2] = remap_field(p[2])
    if len(p) > 3:
        p[3] = remap_field(p[3])
    return "\t".join(p)


def row_faa_rewrite(line: str) -> str:
    p = line.split("\t")
    if len(p) > 1:
        p[1] = remap_header(p[1])
    if len(p) > 2:
        p[2] = remap_header(p[2])
    if len(p) > 4:
        p[4] = remap_token(p[4])
    if len(p) > 5:
        p[5] = remap_token(p[5])
    return "\t".join(p)


# --- io helpers --------------------------------------------------------------
def atomic_write(path: Path, text: str):
    tmp = path.with_name(path.name + TMP)
    tmp.write_text(text)
    os.replace(tmp, path)  # atomic on same filesystem


def backup_once(path: Path):
    bak = path.with_name(path.name + BAK)
    if not bak.exists():
        bak.write_bytes(path.read_bytes())


def existing_baks():
    baks = list(ROOT.glob("*" + BAK))
    baks += list((ROOT / "database" / "faa").glob("*" + BAK))
    return baks


# --- tsv driver (header-aware, terminator-preserving) ------------------------
def process_tsv(path: Path, row_fn, has_header: bool):
    if not path.exists():
        print(f"  SKIP (missing): {path.name}")
        return
    raw = path.read_bytes()
    term = "\r\n" if b"\r\n" in raw else "\n"
    text = raw.decode()
    had_trailing = text.endswith(term)
    body = text[: -len(term)] if had_trailing else text
    lines = body.split(term)
    out, changed = [], 0
    for i, line in enumerate(lines):
        if has_header and i == 0:
            out.append(line)
            continue
        new = row_fn(line)
        if new != line:
            changed += 1
        out.append(new)
    term_label = "CRLF" if term == "\r\n" else "LF"
    print(f"  {path.name}: {changed} rows changed (of {len(lines)}) [{term_label}]")
    if APPLY and changed:
        backup_once(path)
        atomic_write(path, term.join(out) + (term if had_trailing else ""))


# --- faa driver (only rewrite files containing a target; endings preserved) --
def process_faa(faa_dir: Path):
    files = sorted(faa_dir.glob("*.faa"))
    n_files, n_headers = 0, 0
    for fp in files:
        data = fp.read_bytes()
        # prefilter matches the rewrite predicate (PLV__/VP__/PLV|/VP|/bare).
        if b">PLV" not in data and b">VP" not in data:
            continue
        out, changed = [], 0
        for line in data.decode().splitlines(keepends=True):
            if line.startswith(">"):
                stripped = line.rstrip("\r\n")
                ending = line[len(stripped):]
                new_body = remap_header(stripped[1:])
                if ">" + new_body != stripped:
                    changed += 1
                out.append(">" + new_body + ending)
            else:
                out.append(line)
        if changed:
            n_files += 1
            n_headers += changed
            if APPLY:
                backup_once(fp)
                atomic_write(fp, "".join(out))
    print(f"  faa: {n_files} files, {n_headers} headers changed")


def main():
    mode = "APPLY" if APPLY else "DRY-RUN"
    print(f"=== PPV rename ({mode}) on {ROOT} ===")
    if APPLY and not FORCE:
        stale = existing_baks()
        if stale:
            print(f"REFUSING: {len(stale)} existing {BAK} file(s) found (stale/partial prior run).")
            print(f"  e.g. {stale[0].name} ... Remove/restore them, or re-run with --force.")
            sys.exit(2)
    print("TSV files:")
    process_tsv(ROOT / "labels.tsv", row_labels, has_header=False)
    process_tsv(ROOT / "inactive_labels.tsv", row_inactive, has_header=True)
    process_tsv(ROOT / "aliases.tsv", row_aliases, has_header=True)
    process_tsv(ROOT / "label_context.tsv", row_label_context, has_header=True)
    process_tsv(ROOT / "faa_rewrite_map.tsv", row_faa_rewrite, has_header=True)
    print("FAA headers:")
    process_faa(ROOT / "database" / "faa")
    if APPLY:
        print(f"APPLIED. Backups: *{BAK}")
    else:
        print("DRY-RUN. Re-run with --apply to write (atomic; backs up each modified file).")


if __name__ == "__main__":
    main()
