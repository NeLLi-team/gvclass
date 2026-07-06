#!/usr/bin/env python3
"""Build data-driven marker groups from reference-protein overlap.

A GROUP is a set of marker models whose pre-extracted reference FAAs
(resources/database/faa/<model>.faa) share reference proteins, i.e. they are the
same gene family detected by several HMMs. Grouped models are later combined into one
{group}.faa -> one alignment -> one tree -> one nearest-neighbour vote, which removes the
per-model vote over-counting that otherwise fragments classification (e.g. a genome whose
MCP cross-hits 8 PLV profiles casting 8 votes).

Method
------
1. For each model FAA, take the set of reference protein headers (>{label}|{protein}).
2. overlap_coef(A,B) = |refs(A) & refs(B)| / min(|refs(A)|, |refs(B)|)  (robust to size skew).
3. Edge if overlap_coef >= THRESH and shared >= MIN_SHARED. Connected components = groups.
   Conservative THRESH avoids merging distinct families; the worst case of under-grouping is
   the status-quo (separate votes), never a wrong merged tree.
4. Models with no strong family edge stay SINGLETONS (used as-is, one model = one tree).

Output: resources/marker_groups.tsv  (group_name, n_members, total_refs,
annotation_coverage, coherence, dom_annotation, flag, members).

Reproduce:  python scripts/build_marker_groups.py
"""
import os
from collections import defaultdict, Counter

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FAA = os.path.join(REPO, "resources", "database", "faa")
ANN = os.path.join(REPO, "resources", "markers", "annotations.tsv")
OUT = os.path.join(REPO, "src", "config", "marker_groups.tsv")

THRESH = 0.50      # overlap-coefficient threshold for a same-family edge
MIN_SHARED = 30    # minimum shared reference proteins for an edge

# curated names for known viral families (signature member -> group name); unique names
SIG = [
    ("gamadvirusMCP", "mcp_ncldv"), ("GVOGm0003", "mcp_ncldv"),
    ("PLV_MCP_1", "mcp_plv"), ("VP_MCP_1", "mcp_vp"), ("VP_Penton_1", "penton_vp"),
    ("GVOGm0054", "polb"), ("VP_ATPase_1", "atpase_vp"), ("ATPase", "atpase_ncldv_mrya"),
    ("GVOGm0760", "atpase_ncldv_mrya"), ("VLTF3", "vltf3"), ("VLTF2", "vltf2"),
    ("GVOGm0172", "tfiib_cyclin"), ("VP_PRO_1", "protease_vp"),
    ("Mirus_Terminase_merged", "terminase"), ("Mirus_MCP", "mcp_mirus"),
]


def header_set(path):
    # Identity = the full first-token header ">PREFIX__genome|protein" (globally unique).
    # NOT the post-"|" suffix: some headers carry a non-unique numeric protein id
    # (e.g. BAC__g1|2 and BAC__g2|2) that would collide and inflate overlap.
    s = set()
    with open(path, "rb") as fh:
        for line in fh:
            if line[:1] == b">":
                s.add(line[1:].split()[0])
    return s


def main():
    # Exclude previously-generated {group}.faa files (named after group_name in an
    # existing lookup) so grouping is computed over MODEL reference FAAs only.
    group_faa = set()
    if os.path.exists(OUT):
        with open(OUT) as fh:
            for line in fh:
                if line.startswith("#") or line.startswith("group_name\t"):
                    continue
                group_faa.add(line.split("\t", 1)[0] + ".faa")
    files = sorted(f for f in os.listdir(FAA) if f.endswith(".faa") and f not in group_faa)
    models = [f[:-4] for f in files]
    ann = {}
    with open(ANN) as fh:
        next(fh, None)
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                ann[p[0]] = p[1][:46]

    # inverted index header-hash -> models, to count pairwise shared refs (bounded pairs)
    inv = {}
    sizes = {}
    for f, m in zip(files, models):
        s = {hash(h) for h in header_set(os.path.join(FAA, f))}
        sizes[m] = len(s)
        for h in s:
            v = inv.get(h)
            if v is None:
                inv[h] = m
            elif isinstance(v, str):
                inv[h] = [v, m]
            else:
                v.append(m)
    pair = defaultdict(int)
    for v in inv.values():
        if isinstance(v, list) and len(v) <= 40:
            for i in range(len(v)):
                for j in range(i + 1, len(v)):
                    a, b = v[i], v[j]
                    pair[(a, b) if a < b else (b, a)] += 1

    # union-find at THRESH
    parent = {m: m for m in models}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    n_edges = 0
    for (a, b), shared in pair.items():
        if shared < MIN_SHARED:
            continue
        if shared / min(sizes[a], sizes[b]) >= THRESH:
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
            n_edges += 1

    comp = defaultdict(list)
    for m in models:
        comp[find(m)].append(m)
    groups = [sorted(ms, key=lambda m: -sizes[m]) for ms in comp.values() if len(ms) > 1]

    # ---- expert curation from the adversarial group audit (2026-06-21) ----
    # Conservative clustering is imperfect; these are the corrections the audit found:
    # spurious bridge members to drop, cellular over-merges to split, one group to
    # dissolve, and REVIEW_mixed groups confirmed to be true single families.
    DROP_MEMBERS = {"OG5938", "OG7123", "OG5776", "OG1864", "OG27147", "OG16493"}
    DISSOLVE = {"OG5478"}  # a member of a group to dissolve (all members -> singletons)
    SPLITS = {  # seed member -> sub-family member lists (members not listed -> singletons)
        "OG1357": [["OG1357", "OG130", "OG14872"], ["OG276", "OG4821", "OG787", "OG1645", "OG5907"]],
        "OG503": [["OG503", "OG1590", "OG1164", "OG3263"], ["OG8381", "OG527", "OG2605"]],
    }
    FORCE_COMBINE = {"GVOGm0172", "COG0086", "OG1190", "OG1419", "OG203", "OG475"}  # audit-cleared families

    curated, forced = [], set()
    for g in groups:
        gs = set(g)
        if gs & DISSOLVE:
            continue
        seed = next((s for s in SPLITS if s in gs), None)
        if seed:
            for sub in SPLITS[seed]:
                sub = [m for m in sub if m in gs]
                if len(sub) >= 2:
                    curated.append(sorted(sub, key=lambda m: -sizes[m]))
            continue
        g2 = [m for m in g if m not in DROP_MEMBERS]
        if len(g2) >= 2:
            curated.append(g2)
            if gs & FORCE_COMBINE:
                forced.add(tuple(g2))
    groups = sorted(curated, key=lambda g: -sum(sizes[m] for m in g))
    singletons = len(models) - sum(len(g) for g in groups)

    used = set()

    def name_of(members):
        mset = set(members)
        for sig, nm in SIG:
            if sig in mset:
                base = nm
                break
        else:
            base = f"grp_{members[0]}"
        nm = base
        k = 2
        while nm in used:
            nm = f"{base}_{k}"
            k += 1
        used.add(nm)
        return nm

    rows = []
    for g in groups:
        anns = [ann[m] for m in g if ann.get(m)]
        cov = len(anns) / len(g)
        if anns:
            dom, c = Counter(anns).most_common(1)[0]
            coh = c / len(anns)
        else:
            dom, coh = "", 0.0
        # flag only a genuine conflict among >=3 annotated members; audit-cleared
        # families (FORCE_COMBINE) are "ok" despite annotation-label disagreement
        if tuple(g) in forced:
            flag = "ok"
        else:
            flag = "REVIEW_mixed" if len(anns) >= 3 and coh < 0.6 else ("no_annotation" if not anns else "ok")
        rows.append((name_of(g), len(g), sum(sizes[m] for m in g), cov, coh, dom, flag, ";".join(g)))

    with open(OUT, "w") as o:
        o.write("# marker groups (same-family models sharing reference proteins) — built by scripts/build_marker_groups.py\n")
        o.write(f"# overlap_coef>={THRESH}, shared>={MIN_SHARED}; {len(rows)} groups, {sum(r[1] for r in rows)} grouped models, {singletons} singletons (singletons used as-is)\n")
        o.write("# flag: ok | no_annotation (clean but unannotated) | REVIEW_mixed (>=3 annotated members disagree -> check for over-merge)\n")
        o.write("# curation: an adversarial 4-agent group audit (2026-06-21) confirmed every combinable group is a true single\n")
        o.write("# gene family; CURATION in this builder applies its corrections (drop spurious bridge members e.g. OG5938 from\n")
        o.write("# mcp_ncldv; split over-merges grp_OG1357/grp_OG503; dissolve grp_OG5478; force-combine audit-cleared families).\n")
        o.write("# caps models (plv_mcp_caps_*, gamadvirusmcp_caps_*) are not yet in resources/database/faa; once their FAAs are\n")
        o.write("# built from the combined proteome they join (validated 2026-06-21, caps-HMM vs gvclass refs): PLV-class -> mcp_plv,\n")
        o.write("# NCV-like/Mriyavirus -> mcp_ncldv; ~15 novel caps groups (Alpenseevirus/Zephyrvirus/VC124-134/Pam2/Metamonada/...) add new refs.\n")
        o.write("group_name\tn_members\ttotal_refs\tannotation_coverage\tcoherence\tdom_annotation\tflag\tmembers\n")
        for r in rows:
            o.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.2f}\t{r[4]:.2f}\t{r[5]}\t{r[6]}\t{r[7]}\n")

    print(f"edges={n_edges}  models={len(models)}  groups={len(rows)}  grouped={sum(r[1] for r in rows)}  singletons={singletons}")
    print(f"wrote {OUT}")
    print("REVIEW_mixed groups:", [r[0] for r in rows if r[6] == "REVIEW_mixed"] or "none")


if __name__ == "__main__":
    main()
