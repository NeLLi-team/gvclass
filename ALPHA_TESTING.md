# Alpha testing the `gvclass-dev` branch

This branch carries the v2.0 candidate: the capscan/PPV taxonomy update, a
re-curated PPV reference set, a streamlined summary table (see CHANGELOG for the
column changes), **VeryFastTree as the default tree builder**
(IQ-TREE `Q.pfam+R10+F` is opt-in via `--tree-method iqtree`), and the opt-in
`--species-tree` feature. It is wired to a separate **alpha-test resource bundle**
so testers get the updated reference models without waiting for the Zenodo release.
`main` is unchanged and still uses the published Zenodo asset.

## 1. Get the branch

**Fresh clone** — checks out `gvclass-dev` directly via `-b`:

```bash
git clone -b gvclass-dev https://github.com/NeLLi-team/gvclass.git
cd gvclass
curl -fsSL https://pixi.sh/install.sh | bash   # if you don't have pixi
```

**Already have gvclass cloned?** A plain `git pull` leaves you on `main` (the v1
release) — you must switch to the dev branch:

```bash
cd gvclass
git fetch origin gvclass-dev
git checkout gvclass-dev
git pull
```

Either way, **confirm you are on the dev branch before continuing** — this must
print `gvclass-dev`:

```bash
git rev-parse --abbrev-ref HEAD
```

## 2. Resources (auto-downloaded)

On first run (or `pixi run setup-db`) gvclass downloads and verifies the pinned
bundle into `./resources/`:

- URL: `https://dl.newlineages.com/gvclass/resources_v1_7_1.tar.gz`
- Version: `v1.7.1` (capscan/PPV-unified models + protist-A32 PPV-reference re-curation)
- Size: ~2.0 GB &nbsp;·&nbsp; SHA-256 `6f8ca4e0f61e094a7d05669e4024e07db9e3c1813fc07172e25113d362512c14`

The URL/version/sha256 are pinned in `config/gvclass_config.yaml`; the download is
checksum-verified before it is installed. Nothing else to configure.

> gvclass now auto-updates the database when the installed `DB_VERSION` is older
> than the version pinned in the config, so an existing `resources/` (v1.6.0 or
> v1.7.0) is refreshed to v1.7.1 on the next run (interactive: prompted;
> non-interactive: automatic). To force a clean re-download, point gvclass at a
> fresh location:
> `export GVCLASS_DB=/path/to/empty/dir` (or move your old `resources/` aside),
> then run `pixi run setup-db`.

Advanced override (e.g. a local mirror): set `GVCLASS_DB_URL` and
`GVCLASS_DB_SHA256` — a custom URL requires the checksum (or
`GVCLASS_DB_ALLOW_UNVERIFIED=1`).

## 3. What to exercise

```bash
# Standard classification (should match main for non-NCLDV/PPV/MIRUS genomes)
pixi run gvclass example -o example_results --threads 8

# New: per-query species tree + placement taxonomy (default)
pixi run gvclass example -o st_results --species-tree --threads 8

# Optional: also build one combined tree over all queries
pixi run gvclass example -o st_results --species-tree-combined --threads 8

# Tree engine note: VeryFastTree is the default (fast). For higher-accuracy
# (much slower) IQ-TREE trees (Q.pfam+R10+F --fast):
pixi run gvclass example -o iqtree_results --tree-method iqtree --threads 8

# Bootstrap-supported species tree (IQ-TREE ultrafast bootstrap; writes a .contree, slowest):
pixi run gvclass example -o ufb_results --species-tree --tree-method iqtree --iqtree-mode ufboot --threads 8
```

Per NCLDV/PPV/MIRUS genome, `--species-tree` writes `out/species_tree/<query>/`
(`<query>.treefile`, `<query>.partitions.txt`, `species_tree_taxonomy.tsv`) and
adds `species_tree_nn_taxonomy`/`_nn_genome`/`_nn_distance`/`_clade_id` to the
summary. See the README "Species tree" section for details and `--species-tree-trim`.

## 4. Reporting

File issues at https://github.com/NeLLi-team/gvclass/issues — please include the
gvclass branch (`git rev-parse --abbrev-ref HEAD`, should be `gvclass-dev`) and
commit (`git rev-parse --short HEAD`), the resource version
(`cat resources/DB_VERSION`), and the command you ran.
