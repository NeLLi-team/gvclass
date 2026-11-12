# Configuration file

The shipped configuration is `config/gvclass_config.yaml`. CLI flags override YAML values; omitted keys retain their built-in defaults.

## File lookup

`--config` selects a file. A path that does not exist is an error (exit status 2); GVClass does not fall back to another file. Without `--config`, GVClass reads the first of these files that exists:

1. `gvclass_config.yaml` in the working directory.
2. `config/gvclass_config.yaml` in the repository.

Built-in defaults are used if neither exists.

A section (`database`, `pipeline`, `quality`) that is empty or not a mapping is an error that names the section. A missing PyYAML is an error, not a silent switch to the defaults. `threads` must be an integer of at least 1.

## Shipped configuration

```yaml
database:
  path: resources
  cache_path:
  download_url:
  download_version: v2.1.0
  download_sha256:

pipeline:
  tree_method: veryfasttree
  iqtree_mode: fast
  mode_fast: true
  completeness_mode: novelty-aware
  sensitive_mode: true
  contig_screen: true
  contigs_min_length: 10000
  threads: 4

quality:
  min_length: 20000
```

## Database keys

| Key | Default | Meaning |
| --- | --- | --- |
| `path` | `resources` | Database directory. Relative paths resolve from the repository root. |
| `cache_path` | unset | Cache for TSV/FASTA files materialized from Parquet. Relative paths resolve inside the database directory. |
| `download_url` | unset | URL of a single resource archive, used with `download_sha256` as a mirror. When unset, setup uses the built-in source: the v2.1.0 archive in 24 parts on Zenodo, each part and the joined archive checked against pinned SHA-256 values. |
| `download_version` | `v2.1.0` | Pinned database version. An older installed `DB_VERSION` produces a prompt (terminal) or a warning (no terminal); `setup-db --update` installs it. |
| `download_sha256` | unset | Expected checksum of the archive at `download_url`. Required when `download_url` is set. |

Database path precedence is `--database`, then `GVCLASS_DB`, then `database.path`.

Setup uses one download source, chosen in this order: the `GVCLASS_DB_URL` environment variable, the `database.download_*` keys, then the built-in v2.1.0 source. An environment URL requires `GVCLASS_DB_SHA256`. The optional `GVCLASS_DB_VERSION` value describes the same archive. Setup does not follow Zenodo's latest link or fall back to an older bundle. A Zenodo archive that GVClass finds automatically must have a valid SHA-256 sidecar. Set `GVCLASS_DB_ALLOW_UNVERIFIED=1` to opt in to an unverified environment archive.

Compact resource bundles can store labels and reference proteins as Parquet (`parquet/labels/labels.parquet` and `parquet/faa.parquet`). Validation accepts reference proteins from `parquet/faa.parquet` or from nonempty legacy `database/faa/*.faa` files.

`GVCLASS_RESOURCE_CACHE` overrides `cache_path`. The default is `<database.path>/.gvclass_cache`. When that directory is unwritable, GVClass uses `<system temp>/gvclass_resource_cache_<uid>/<database-path hash>`, which it creates with mode 0700 and refuses to use if another user owns it. The cache keeps one subdirectory per resource identity (the resource manifest digest, or the resolved database path when there is no manifest, and the size and modification time of each Parquet table), so a replaced resource never reads files cached from an older one, and two resources can share one cache directory. Remove old subdirectories to free space. See [configure the database](../how-to/configure-the-database.md).

## Pipeline keys

| Key | Default | Meaning |
| --- | --- | --- |
| `tree_method` | `veryfasttree` | Tree builder: `veryfasttree`, `iqtree`, or `fasttree` (alias for `veryfasttree`). |
| `iqtree_mode` | `fast` | Species-tree IQ-TREE search: `fast` or `ufboot`. Per-marker trees use `--fast`. |
| `mode_fast` | `true` | Skip order-level marker trees; retain order-level HMM searches. |
| `completeness_mode` | `novelty-aware` | Primary completeness estimator: `novelty-aware` or `legacy`. |
| `sensitive_mode` | `true` | Use `E=1e-5` and `domE=1e-5` instead of GA cutoffs. |
| `contig_screen` | `true` | Screen each nucleotide bin after the queries finish, when the add-on is installed. A bin that the cellular screen scores takes the cellular share as `estimated_contamination`. `--no-screen` turns it off for one run. |
| `contigs_min_length` | `10000` | Minimum contig length in bp for `--contigs`. |
| `threads` | `4` | Total thread budget, an integer of at least 1. |

## Quality keys

| Key | Default | Meaning |
| --- | --- | --- |
| `min_length` | `20000` | Minimum total nucleotide length in bp per input in bin/MAG mode. Overridden by `--min-length` or bypassed by `--allow-short`. |

Genetic codes are fixed: `0, 1, 4, 6, 11, 15, 29, 106, 129`. See [markers and genetic codes](markers.md). Log verbosity is set with `-v`/`--verbose`.
