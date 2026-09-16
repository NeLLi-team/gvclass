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
  download_url: https://zenodo.org/records/21225457/files/resources_v2_0_0.tar.gz?download=1
  download_version: v2.0.0
  download_sha256: df1c3a9d15a90307775f42f57e2a7c89436ed523883025f6fc94013035f5e066

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
| `download_version` | `v2.0.0` | Pinned database version. An older installed `DB_VERSION` triggers a download. |
| `download_sha256` | SHA-256 above | Expected archive checksum. |

Database path precedence is `--database`, then `GVCLASS_DB`, then `database.path`.

`GVCLASS_RESOURCE_CACHE` overrides `cache_path`. The default is `<database.path>/.gvclass_cache`, with a hashed system-temporary cache used when the database directory is unwritable. See [configure the database](../how-to/configure-the-database.md).

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
