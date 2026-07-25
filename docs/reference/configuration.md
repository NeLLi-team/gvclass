# Configuration file

GVClass reads run settings from a YAML file. Command-line flags override file values, and file values override built-in defaults. The shipped file is `config/gvclass_config.yaml`.

## File lookup and precedence

GVClass resolves the config in this order:

1. The path passed to `-c` / `--config`.
2. `./gvclass_config.yaml` in the current working directory.
3. `config/gvclass_config.yaml` in the repository.

Effective values are merged with this precedence:

| Source | Precedence |
| --- | --- |
| CLI flags (see [cli.md](cli.md)) | highest |
| Config file | middle |
| Built-in defaults | lowest |

## Shipped configuration

The block below is the public setup-download pin for the v2.0.0 resource bundle.

```yaml
# Database configuration
database:
  path: resources
  cache_path:
  download_url: https://zenodo.org/records/21225457/files/resources_v2_0_0.tar.gz?download=1
  download_version: v2.0.0
  download_sha256: df1c3a9d15a90307775f42f57e2a7c89436ed523883025f6fc94013035f5e066

# Pipeline settings
pipeline:
  tree_method: veryfasttree
  iqtree_mode: fast
  mode_fast: true
  completeness_mode: novelty-aware
  sensitive_mode: true
  contigs_min_length: 10000
  threads: 4
  output_pattern: "{query_dir}_results"

# Quality thresholds
quality:
  min_length: 20000
```

Every key GVClass reads is listed above. The file is the single source of the
built-in defaults: `src/bin/gvclass_cli.py` carries a matching literal for
installs where `config/` is absent, and a test fails if the two ever disagree.

## Keys

### database

| Key | Default | Meaning |
| --- | --- | --- |
| `path` | `resources` | Database bundle location. Relative paths resolve from the repo root. |
| `cache_path` | unset | Optional materialization cache for compact Parquet resources. Relative paths resolve inside the selected database directory. `GVCLASS_RESOURCE_CACHE` overrides this key. When unset, GVClass uses `<database.path>/.gvclass_cache` and falls back to a hashed system-temp cache only if the database directory is not writable. |
| `download_url` | `https://zenodo.org/records/21225457/files/resources_v2_0_0.tar.gz?download=1` | Archive fetched for first-time setup and auto-update. |
| `download_version` | `v2.0.0` | Pinned public download version. An older installed `DB_VERSION` triggers a re-download. |
| `download_sha256` | `df1c3a9d...35f5e066` | Checksum verified after download. |

The full checksum is `df1c3a9d15a90307775f42f57e2a7c89436ed523883025f6fc94013035f5e066`. Override the database location with `-d` / `--database` or the `GVCLASS_DB` environment variable. See [configure the database](../how-to/configure-the-database.md).

Compact resource bundles can store labels and reference proteins as Parquet (`parquet/labels/labels.parquet` and `parquet/faa.parquet`) instead of `labels.tsv` and `database/faa/*.faa`. GVClass materializes the legacy views it needs into `<database.path>/.gvclass_cache/` by default. Set `database.cache_path` or `GVCLASS_RESOURCE_CACHE` to move that cache outside the database directory; the environment variable has precedence.

### pipeline

| Key | Default | Meaning |
| --- | --- | --- |
| `tree_method` | `veryfasttree` | Per-marker tree builder. One of `veryfasttree`, `iqtree`, `fasttree` (`fasttree` is an alias for `veryfasttree`). |
| `iqtree_mode` | `fast` | Species-tree IQ-TREE search. `fast` or `ufboot`. Per-marker trees always run `--fast`. |
| `mode_fast` | `true` | Skip order-level marker trees when `true`. |
| `completeness_mode` | `novelty-aware` | Completeness estimator surfaced in the summary. `novelty-aware` or `legacy`. |
| `sensitive_mode` | `true` | Use `E=1e-5` cutoffs instead of GA model cutoffs. |
| `contigs_min_length` | `10000` | Minimum contig length (bp) when splitting files in `--contigs` mode. Override per run with `--contigs-min-length`. |
| `threads` | `4` | Total compute threads. |
| `output_pattern` | `{query_dir}_results` | Output directory name. `{query_dir}` is the basename of the query directory. |

For tuning these against runtime, see [tune speed and accuracy](../how-to/tune-speed-and-accuracy.md).

### quality

| Key | Default | Meaning |
| --- | --- | --- |
| `min_length` | `20000` | Minimum total nucleotide length (bp) for each `.fna` input in bin/MAG mode. Override per run with `--min-length`; use `--allow-short` to bypass the gate. |

## Genetic codes are not configurable

GVClass always evaluates the full panel `0, 1, 4, 6, 11, 15, 29, 106, 129` for
every nucleotide query and keeps the best-scoring code. Code `0` is pyrodigal
meta mode. See [markers and genetic codes](markers.md).

## Log verbosity

Log level follows `-v` / `--verbose` on the command line, not a config key.
