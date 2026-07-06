# Configure the database

GVClass reads its reference markers, HMMs, and tree references from one database directory. You control where that directory lives, point several runs at a single shared copy, and let GVClass keep it current. The bundle is versioned separately from the software. The public setup-download block currently installs `DB_VERSION` v2.0.0.

## Where GVClass looks for the database

GVClass resolves the database location from four sources. The first one that is set wins.

| Priority | Source | How to set |
|----------|--------|------------|
| 1 | `--database` flag (`-d`) | `pixi run gvclass <in> -o <out> -d /path/to/gvclass_resources` |
| 2 | `GVCLASS_DB` environment variable | `export GVCLASS_DB=/path/to/gvclass_resources` |
| 3 | `database.path` in the config | edit `config/gvclass_config.yaml` |
| 4 | default `<repo>/resources` | nothing; used when none of the above is set |

Relative paths resolve against the repo directory, so the default `path: resources` becomes `<repo>/resources`. See the [CLI reference](../reference/cli.md) for the `--database` flag and the [configuration reference](../reference/configuration.md) for every key.

## Set up the database the first time

Run the setup task once per database location. It downloads the bundle, verifies the SHA-256 and size, extracts it, and writes a `DB_VERSION` file.

```bash
pixi run setup-db
```

The download is about 1.5 GB. Once it finishes, a normal run finds the database under `<repo>/resources` with no extra flags. For a full first run, follow the [getting-started tutorial](../tutorials/getting-started.md).

## Share one database across runs

On shared storage, download the bundle once and point every run at it. Set `GVCLASS_DB` to the shared path, run setup there, then run as usual.

```bash
export GVCLASS_DB=/path/to/shared/gvclass_resources
pixi run setup-db
pixi run gvclass <input_dir> -o <output_dir> -t 16
```

Every shell that exports the same `GVCLASS_DB` reuses that copy, which avoids a per-user 1.5 GB download. This is the pattern to use on an HPC scratch or project filesystem; see [Run on HPC](run-on-hpc.md).

!!! tip
    Put the `export GVCLASS_DB=...` line in your shell profile or job script so every run picks it up automatically.

## The database config block

The `database` block of `config/gvclass_config.yaml` controls the path and the bundle that setup downloads. These values describe the public download archive, not necessarily the already-installed repo-local `resources/` tree.

```yaml
database:
  path: resources
  cache_path:
  download_url: https://zenodo.org/records/21225457/files/resources_v2_0_0.tar.gz?download=1
  download_version: v2.0.0
  download_sha256: df1c3a9d15a90307775f42f57e2a7c89436ed523883025f6fc94013035f5e066
  expected_size: 1497
```

Each key is defined in the [configuration reference](../reference/configuration.md); the values above are the current public setup-download pin. The keys you are most likely to change are `path` (where the database lives), `cache_path` (where compact-resource views are materialized), and `download_version` (the pinned download bundle, currently `v2.0.0`).

The v2.0.0 archive is published on Zenodo as DOI [10.5281/zenodo.21225457](https://doi.org/10.5281/zenodo.21225457).

## Label namespaces in the database

The selected database bundle controls the reference labels that GVClass can emit. The v2.0.0 repo-local resources include putative endogenous viral element references as `EUK-pEVE__...` FASTA and label IDs, with taxonomy strings that carry `-pEVE` at every rank. No command-line flag is required to activate that namespace; point GVClass at the bundle with `--database`, `GVCLASS_DB`, or the `database.path` config key. In species-tree runs, `EUK-pEVE__...` references are eligible auxiliary leaves for the NCLDV, PPV, and MIRUS panels, but ordinary `EUK__...` references are not.

For interpretation of `EUK-pEVE` outputs, see [Taxonomy and classification](../explanation/taxonomy.md#putative-eve-references).

## Compact Parquet resources

GVClass supports two on-disk layouts for labels and reference proteins:

- Legacy bundles store active labels in `labels.tsv` and marker reference proteins in `database/faa/<marker>.faa`.
- Compact bundles may instead store active labels in `parquet/labels/labels.parquet` and all marker proteins in `parquet/faa.parquet`.

The compact EUK80 keep-pEVE projection collapses ordinary `EUK__...` reference proteins to 80% amino-acid identity representatives while keeping all `EUK-pEVE__...` proteins. The taxonomy namespace does not change: ordinary eukaryotic labels remain `EUK__...`, and pEVE labels remain `EUK-pEVE__...`.

When a compact bundle is used, GVClass materializes only the label table and marker FASTA files needed for the current run into `.gvclass_cache/` inside the database directory. The archive does not ship cache contents; a copied resource directory warms its cache beside the bundle on first use. Set `database.cache_path` or `GVCLASS_RESOURCE_CACHE=/path/to/cache` if the database directory is read-only or if you want the materialized files elsewhere; the environment variable takes precedence over the config key.

## Keep the database current

GVClass compares the installed `DB_VERSION` against the configured `download_version` at the start of a run. When the installed bundle is older, it re-downloads:

- Interactive runs prompt before updating, with yes as the default answer.
- Non-interactive runs (scripts, Slurm jobs) update automatically rather than run against a stale bundle.

This works for any configured source, including the non-Zenodo Cloudflare bundle. To force a fresh download at any time, delete the database directory and run `pixi run setup-db` again.

For the full set of configuration keys, see the [configuration reference](../reference/configuration.md).
