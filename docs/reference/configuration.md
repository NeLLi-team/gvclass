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

## Shipped configuration (gvclass-dev)

The block below is the public setup-download pin. The checked repo-local
`resources/` tree can be newer; current local resources are `DB_VERSION`
v2.0.0, while this download block stays on v1.7.1 until a matching v2.0.0
archive and checksum are published.

```yaml
# Database configuration
database:
  path: resources
  cache_path:
  download_url: https://dl.newlineages.com/gvclass/resources_v1_7_1.tar.gz
  download_version: v1.7.1
  download_sha256: 6f8ca4e0f61e094a7d05669e4024e07db9e3c1813fc07172e25113d362512c14
  expected_size: 2005

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

# Genetic code settings
genetic_codes:
  codes: [0, 1, 4, 6, 11, 15, 29, 106, 129]
  improvement_threshold: 0.05

# Quality thresholds
quality:
  min_markers: 3
  min_length: 20000
  recommended_length: 50000

# Resource allocation
resources:
  memory_limit: 8
  task_timeout: 60

# Logging
logging:
  level: ERROR
  keep_temp: false
```

!!! note
    On `main`, the `database` block points at the Zenodo release asset instead of the `dl.newlineages.com` tunnel shown here. See [configure the database](../how-to/configure-the-database.md).

## Keys

### database

| Key | Default | Meaning |
| --- | --- | --- |
| `path` | `resources` | Database bundle location. Relative paths resolve from the repo root. |
| `cache_path` | unset | Optional materialization cache for compact Parquet resources. Relative paths resolve inside the selected database directory. `GVCLASS_RESOURCE_CACHE` overrides this key. When unset, GVClass uses `<database.path>/.gvclass_cache` and falls back to a hashed system-temp cache only if the database directory is not writable. |
| `download_url` | `https://dl.newlineages.com/gvclass/resources_v1_7_1.tar.gz` | Archive fetched for first-time setup and auto-update; intentionally left at v1.7.1 until the v2.0.0 archive exists. |
| `download_version` | `v1.7.1` | Pinned public download version. An older installed `DB_VERSION` triggers a re-download. |
| `download_sha256` | `6f8ca4e0...362512c14` | Checksum verified after download. |
| `expected_size` | `2005` | Expected archive size in MB, used for validation. |

The full checksum is `6f8ca4e0f61e094a7d05669e4024e07db9e3c1813fc07172e25113d362512c14`. Override the database location with `-d` / `--database` or the `GVCLASS_DB` environment variable. See [configure the database](../how-to/configure-the-database.md).

Compact resource bundles can store labels and reference proteins as Parquet (`parquet/labels/labels.parquet` and `parquet/faa.parquet`) instead of `labels.tsv` and `database/faa/*.faa`. GVClass materializes the legacy views it needs into `<database.path>/.gvclass_cache/` by default. Set `database.cache_path` or `GVCLASS_RESOURCE_CACHE` to move that cache outside the database directory; the environment variable has precedence.

### pipeline

| Key | Default | Meaning |
| --- | --- | --- |
| `tree_method` | `veryfasttree` | Per-marker tree builder. One of `veryfasttree`, `iqtree`, `fasttree` (`fasttree` is an alias for `veryfasttree`). |
| `iqtree_mode` | `fast` | Species-tree IQ-TREE search. `fast` or `ufboot`. Per-marker trees always run `--fast`. |
| `mode_fast` | `true` | Skip order-level marker trees when `true`. |
| `completeness_mode` | `novelty-aware` | Completeness estimator surfaced in the summary. `novelty-aware` or `legacy`. |
| `sensitive_mode` | `true` | Use `E=1e-5` cutoffs instead of GA model cutoffs. |
| `contigs_min_length` | `10000` | Minimum contig length (bp) when splitting files in `--contigs` mode. |
| `threads` | `4` | Total compute threads. |
| `output_pattern` | `{query_dir}_results` | Output directory name. `{query_dir}` is the basename of the query directory. |

For tuning these against runtime, see [tune speed and accuracy](../how-to/tune-speed-and-accuracy.md).

### genetic_codes

| Key | Default | Meaning |
| --- | --- | --- |
| `codes` | `[0, 1, 4, 6, 11, 15, 29, 106, 129]` | Genetic codes tested during gene calling. Code `0` is pyrodigal meta mode. |
| `improvement_threshold` | `0.05` | Minimum fractional improvement (5%) required to select an alternative code over meta. |

### quality

| Key | Default | Meaning |
| --- | --- | --- |
| `min_markers` | `3` | Minimum markers required for a query. |
| `min_length` | `20000` | Minimum contig length (bp). Shorter `.fna` inputs are skipped unless `--allow-short` is set. |
| `recommended_length` | `50000` | Recommended minimum length (bp). |

### resources

| Key | Default | Meaning |
| --- | --- | --- |
| `memory_limit` | `8` | Memory limit per worker (GB). |
| `task_timeout` | `60` | Per-task timeout (minutes). |

### logging

| Key | Default | Meaning |
| --- | --- | --- |
| `level` | `ERROR` | Log level. One of `DEBUG`, `INFO`, `WARNING`, `ERROR`, `CRITICAL`. |
| `keep_temp` | `false` | Keep intermediate files when `true`. |
