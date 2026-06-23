#!/usr/bin/env python3
"""Combine each marker group's member reference FAAs into one {group}.faa.

For every combinable group in src/config/marker_groups.tsv (flag ok|no_annotation),
concatenate the member models' resources/database/faa/<model>.faa, de-duplicating
proteins by header so a protein shared by several member models appears ONCE. The
combined {group}.faa is the reference set for the group's single tree (one group ->
one alignment -> one tree -> one consolidated nearest-neighbour vote).

Standalone (no src import, so it runs without the runtime deps).
Usage: python scripts/build_group_faas.py [--groups a,b,c] [--apply]
       (default: dry-run, all combinable groups)
"""
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAA = os.path.join(REPO, "resources", "database", "faa")
LOOKUP = os.path.join(REPO, "src", "config", "marker_groups.tsv")
APPLY = "--apply" in sys.argv

only = None
for i, a in enumerate(sys.argv):
    if a == "--groups" and i + 1 < len(sys.argv):
        only = set(sys.argv[i + 1].split(","))


def load_groups(path, flags=frozenset({"ok", "no_annotation"})):
    g2m = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("group_name\t"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) >= 8 and p[6] in flags:
                members = [m for m in p[7].split(";") if m]
                if len(members) >= 2:
                    g2m[p[0]] = members
    return g2m


def main():
    groups = load_groups(LOOKUP)
    if only:
        groups = {g: m for g, m in groups.items() if g in only}
    print(f"{'group':<22}{'models':<8}{'raw_seqs':<10}{'unique':<10}redundant_removed")
    for group, models in sorted(groups.items(), key=lambda kv: kv[0]):
        seen, recs, raw, nmodel = set(), [], 0, 0
        for m in models:
            p = os.path.join(FAA, f"{m}.faa")
            if not os.path.exists(p):
                print(f"  WARN missing {m}.faa")
                continue
            nmodel += 1
            keep = False
            with open(p) as fh:
                for line in fh:
                    if line.startswith(">"):
                        raw += 1
                        h = line[1:].split()[0]
                        keep = h not in seen
                        if keep:
                            seen.add(h)
                            recs.append(line)
                    elif keep:
                        recs.append(line)
        uniq = len(seen)
        print(f"{group:<22}{nmodel:<8}{raw:<10}{uniq:<10}{raw - uniq}")
        if APPLY:
            with open(os.path.join(FAA, f"{group}.faa"), "w") as o:
                o.writelines(recs)
    if not APPLY:
        print("\nDRY-RUN — add --apply to write resources/database/faa/{group}.faa")


if __name__ == "__main__":
    main()
