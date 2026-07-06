#!/usr/bin/env python3
"""Independent post-apply verification of the PPV rename (NOT using the migration code).

Compares each *.bak_ppv (pre) against the migrated file (post) and asserts:
  * row count unchanged
  * conservation:  PPV__(post) == PLV__(pre) + VP__(pre)   [token-start occurrences]
                   'PPV|'(post) == 'PLV|'(pre) + 'VP|'(pre)
  * lifestyle tags PLV_unclassified / VP_unclassified counts UNCHANGED
  * zero residual domain tokens post: 'PLV|','VP|', and token-start 'PLV__','VP__'
And for faa: every header genome-id resolves in labels.tsv col0 (bijection sanity).

Exit 0 = all pass. Usage: python scripts/verify_ppv_migration.py
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "resources"
BAK = ".bak_ppv"

# token-start occurrences: start of string, or preceded by tab '|' '>' ';'
TOKSTART = r"(?:^|[\t|>;])"
RE_PLV_ID = re.compile(TOKSTART + r"PLV__")
RE_VP_ID = re.compile(TOKSTART + r"VP__")
RE_PPV_ID = re.compile(TOKSTART + r"PPV__")
RE_PLV_DOM = re.compile(TOKSTART + r"PLV\|")
RE_VP_DOM = re.compile(TOKSTART + r"VP\|")
RE_PPV_DOM = re.compile(TOKSTART + r"PPV\|")
RE_PLV_LIFE = re.compile(r"PLV_unclassified")
RE_VP_LIFE = re.compile(r"VP_unclassified")

fails = []


def count(text, rx):
    return len(rx.findall(text))


def check_tsv(name):
    post = ROOT / name
    pre = ROOT / (name + BAK)
    if not post.exists():
        print(f"  {name}: SKIP (missing)")
        return
    if not pre.exists():
        print(f"  {name}: NOTE no backup (unchanged by migration) — checking residuals only")
        t = post.read_text()
        res = count(t, RE_PLV_ID) + count(t, RE_VP_ID) + count(t, RE_PLV_DOM) + count(t, RE_VP_DOM)
        if res:
            fails.append(f"{name}: {res} residual PLV/VP domain tokens with no backup")
        print(f"  {name}: residual PLV/VP domain tokens = {res} (want 0)")
        return
    a, b = pre.read_text(), post.read_text()
    pre_rows, post_rows = a.count("\n"), b.count("\n")
    plv_id, vp_id = count(a, RE_PLV_ID), count(a, RE_VP_ID)
    ppv_id_post = count(b, RE_PPV_ID)
    plv_dom, vp_dom = count(a, RE_PLV_DOM), count(a, RE_VP_DOM)
    ppv_dom_post = count(b, RE_PPV_DOM)
    # residuals post
    res_id = count(b, RE_PLV_ID) + count(b, RE_VP_ID)
    res_dom = count(b, RE_PLV_DOM) + count(b, RE_VP_DOM)
    # lifestyle conservation
    life_ok = count(a, RE_PLV_LIFE) == count(b, RE_PLV_LIFE) and count(a, RE_VP_LIFE) == count(b, RE_VP_LIFE)

    ok = True
    if pre_rows != post_rows:
        ok = False; fails.append(f"{name}: row count {pre_rows}->{post_rows}")
    # conservation (PPV gained == PLV+VP lost), accounting for any pre-existing PPV
    ppv_id_pre = count(a, RE_PPV_ID); ppv_dom_pre = count(a, RE_PPV_DOM)
    if ppv_id_post - ppv_id_pre != plv_id + vp_id:
        ok = False; fails.append(f"{name}: ID conservation off (PPV+{ppv_id_post-ppv_id_pre} vs PLV+VP {plv_id+vp_id})")
    if ppv_dom_post - ppv_dom_pre != plv_dom + vp_dom:
        ok = False; fails.append(f"{name}: domain-token conservation off (PPV+{ppv_dom_post-ppv_dom_pre} vs {plv_dom+vp_dom})")
    if res_id or res_dom:
        ok = False; fails.append(f"{name}: residual PLV/VP tokens post (id={res_id} dom={res_dom})")
    if not life_ok:
        ok = False; fails.append(f"{name}: lifestyle-tag counts changed")
    print(f"  {name}: rows {post_rows} | PLV__+VP__ {plv_id+vp_id}->PPV__ {ppv_id_post-ppv_id_pre} | "
          f"dom {plv_dom+vp_dom}->{ppv_dom_post-ppv_dom_pre} | residual id={res_id} dom={res_dom} | "
          f"lifestyle {'ok' if life_ok else 'CHANGED'} | {'OK' if ok else 'FAIL'}")


def _resolves(gid, label_ids):
    # mirror the runtime resolvers (contamination_scoring._parse_subject_genome_id /
    # summarize_full._extract_genome_id): exact id, else strip a trailing _<int> gene index.
    if gid in label_ids:
        return True
    return re.sub(r"_\d+$", "", gid) in label_ids


def check_faa_bijection():
    labels = ROOT / "labels.tsv"
    label_ids = set()
    for line in labels.read_text().splitlines():
        if line:
            label_ids.add(line.split("\t", 1)[0])
    faa_dir = ROOT / "database" / "faa"
    missing, checked, residual = set(), 0, 0
    for fp in faa_dir.glob("*.faa"):
        data = fp.read_bytes()
        if b">PLV__" in data or b">VP__" in data:
            residual += 1
        if b">PPV__" not in data:
            continue
        for line in data.decode().splitlines():
            if line.startswith(">PPV__"):
                gid = line[1:].split("|", 1)[0]
                checked += 1
                if not _resolves(gid, label_ids):
                    missing.add(gid)
    print(f"  faa: PPV headers checked={checked}, files with residual >PLV__/>VP__={residual}, "
          f"genome-ids missing from labels={len(missing)}")
    if residual:
        fails.append(f"faa: {residual} files still have >PLV__/>VP__ headers")
    if missing:
        fails.append(f"faa: {len(missing)} PPV genome-ids absent from labels.tsv e.g. {list(missing)[:3]}")


def main():
    print("=== PPV migration verification (independent of migration code) ===")
    print("TSV conservation/residual/lifestyle:")
    for n in ["labels.tsv", "inactive_labels.tsv", "aliases.tsv", "label_context.tsv", "faa_rewrite_map.tsv"]:
        check_tsv(n)
    print("FAA bijection:")
    check_faa_bijection()
    print()
    if fails:
        print("*** VERIFICATION FAILED ***")
        for f in fails:
            print("  -", f)
        sys.exit(1)
    print("ALL VERIFICATION CHECKS PASSED")


if __name__ == "__main__":
    main()
