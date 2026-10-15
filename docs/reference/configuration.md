# Configuration file

The shipped configuration is `config/gvclass_config.yaml`. CLI flags override YAML values; omitted keys retain their built-in defaults.

## File lookup

`--config` selects a file; its default is `config/gvclass_config.yaml`. An absolute path is read directly. For a relative path, GVClass checks:

1. The supplied path relative to the working directory.
2. `gvclass_config.yaml` in the working directory.
3. The supplied path relative to the repository root.

Built-in defaults are used if no file is found.

## Shipped configuration

```yaml
database:
  path: resources
  cache_path:
  download_url: https://zenodo.org/records/23381313/files/resources_v2_1_0.tar.gz?download=1
  download_version: v2.1.0
  download_sha256: 823023a0d9af82d6b01d7beeb1a12af422432c7c1e3a11b57f04130bfa7f5ee6

pipeline:
  tree_method: veryfasttree
  iqtree_mode: fast
  mode_fast: true
  completeness_mode: novelty-aware
  sensitive_mode: true
  contigs_min_length: 10000
  threads: 4
  output_pattern: "{query_dir}_results"

quality:
  min_length: 20000
```

## Database keys

| Key | Default | Meaning |
| --- | --- | --- |
| `path` | `resources` | Database directory. Relative paths resolve from the repository root. |
| `cache_path` | unset | Cache for TSV/FASTA files materialized from Parquet. Relative paths resolve inside the database directory. |
| `download_url` | URL above | Resource archive downloaded during setup or an update. |
| `download_version` | `v2.1.0` | Pinned database version. An older installed `DB_VERSION` triggers a download. |
| `download_sha256` | SHA-256 above | Expected archive checksum. |

Database path precedence is `--database`, then `GVCLASS_DB`, then `database.path`.

Setup uses one download source, chosen in this order: the `GVCLASS_DB_URL` environment variable, the `database.download_*` keys, then the built-in v2.1.0 source. An environment URL requires `GVCLASS_DB_SHA256`. The optional `GVCLASS_DB_VERSION` and `GVCLASS_DB_FILENAME` values describe the same archive. Setup uses only the selected source. It does not follow Zenodo's latest link or fall back to an older bundle. A Zenodo archive that GVClass finds automatically must have a valid SHA-256 sidecar. Set `GVCLASS_DB_ALLOW_UNVERIFIED=1` to opt in to an unverified environment archive.

Compact resource bundles can store labels and reference proteins as Parquet (`parquet/labels/labels.parquet` and `parquet/faa.parquet`). Validation accepts reference proteins from `parquet/faa.parquet` or from nonempty legacy `database/faa/*.faa` files.

`GVCLASS_RESOURCE_CACHE` overrides `cache_path`. The default is `<database.path>/.gvclass_cache`, with a hashed system-temporary cache used when the database directory is unwritable. The cache keeps one subdirectory per resource identity (the resource manifest digest and the size and modification time of each Parquet table), so a replaced resource never reads files cached from an older one, and two resources can share one cache directory. Remove old subdirectories to free space. See [configure the database](../how-to/configure-the-database.md).

## Pipeline keys

| Key | Default | Meaning |
| --- | --- | --- |
| `tree_method` | `veryfasttree` | Tree builder: `veryfasttree`, `iqtree`, or `fasttree` (alias for `veryfasttree`). |
| `iqtree_mode` | `fast` | Species-tree IQ-TREE search: `fast` or `ufboot`. Per-marker trees use `--fast`. |
| `mode_fast` | `true` | Skip order-level marker trees; retain order-level HMM searches. |
| `completeness_mode` | `novelty-aware` | Primary completeness estimator: `novelty-aware` or `legacy`. |
| `sensitive_mode` | `true` | Use `E=1e-5` and `domE=1e-5` instead of GA cutoffs. |
| `contigs_min_length` | `10000` | Minimum contig length in bp for `--contigs`. |
| `threads` | `4` | Total thread budget. |
| `output_pattern` | `{query_dir}_results` | Default output directory. `{query_dir}` is the query directory's basename. |

## Quality keys

| Key | Default | Meaning |
| --- | --- | --- |
| `min_length` | `20000` | Minimum total nucleotide length in bp per input in bin/MAG mode. Overridden by `--min-length` or bypassed by `--allow-short`. |

Genetic codes are fixed: `0, 1, 4, 6, 11, 15, 29, 106, 129`. See [markers and genetic codes](markers.md). Log verbosity is set with `-v`/`--verbose`.
