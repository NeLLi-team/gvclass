# Configure the database

GVClass reads its reference markers, HMMs, and tree references from one database directory. You control where that directory lives, point several runs at a single shared copy, and let GVClass keep it current. The bundle is versioned separately from the software; this page covers the v1.7.1 bundle that ships with GVClass 2.0.0.

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

The download is about 2 GB. Once it finishes, a normal run finds the database under `<repo>/resources` with no extra flags. For a full first run, follow the [getting-started tutorial](../tutorials/getting-started.md).

## Share one database across runs

On shared storage, download the bundle once and point every run at it. Set `GVCLASS_DB` to the shared path, run setup there, then run as usual.

```bash
export GVCLASS_DB=/path/to/shared/gvclass_resources
pixi run setup-db
pixi run gvclass <input_dir> -o <output_dir> -t 16
```

Every shell that exports the same `GVCLASS_DB` reuses that copy, which avoids a per-user 2 GB download. This is the pattern to use on an HPC scratch or project filesystem; see [Run on HPC](run-on-hpc.md).

!!! tip
    Put the `export GVCLASS_DB=...` line in your shell profile or job script so every run picks it up automatically.

## The database config block

The `database` block of `config/gvclass_config.yaml` controls the path and the bundle that setup downloads.

```yaml
database:
  path: resources
  download_url: https://dl.newlineages.com/gvclass/resources_v1_7_1.tar.gz
  download_version: v1.7.1
  download_sha256: 6f8ca4e0f61e094a7d05669e4024e07db9e3c1813fc07172e25113d362512c14
  expected_size: 2005
```

Each key is defined in the [configuration reference](../reference/configuration.md); the values above are what ships on the `gvclass-dev` branch. The two you are most likely to change are `path` (where the database lives) and `download_version` (the pinned bundle, currently `v1.7.1`).

!!! note
    On the `gvclass-dev` branch, `download_url` points at the Cloudflare tunnel bundle (`dl.newlineages.com`, v1.7.1). The `main` branch points at the Zenodo asset instead. Both serve the same v1.7.1 contents.

## Keep the database current

GVClass compares the installed `DB_VERSION` against the configured `download_version` at the start of a run. When the installed bundle is older, it re-downloads:

- Interactive runs prompt before updating, with yes as the default answer.
- Non-interactive runs (scripts, Slurm jobs) update automatically rather than run against a stale bundle.

This works for any configured source, including the non-Zenodo Cloudflare bundle. To force a fresh download at any time, delete the database directory and run `pixi run setup-db` again.

For the full set of configuration keys, see the [configuration reference](../reference/configuration.md).
