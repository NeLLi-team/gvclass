# Configure the database

Run these commands from the GVClass repository to install or select a reference
database. Database and software versions are separate.

## Install in the default location

```bash
pixi run setup-db
cat resources/DB_VERSION
```

Setup downloads and verifies the archive, then extracts it into `resources/`.
The configured bundle is about 1.5 GB and is published on
[Zenodo](https://doi.org/10.5281/zenodo.21225457). `DB_VERSION` records the installed
version.

## Use a shared database

1. Set `GVCLASS_DB` to an absolute path. Replace the example path with a directory
   on your shared filesystem.

    ```bash
    export GVCLASS_DB=/path/to/shared/gvclass_resources
    ```

2. Set up the database at that location. If a complete, current copy exists,
   setup reuses it.

    ```bash
    pixi run setup-db
    cat "$GVCLASS_DB/DB_VERSION"
    ```

3. Keep the variable set in each shell or job script that runs GVClass.

    ```bash
    pixi run gvclass example -o example_shared_results -t 8
    ```

To select the database for one command, use
`--database /path/to/shared/gvclass_resources`. Check the selected path in the
startup output.

Database location is resolved in this order:

| Priority | Setting |
| --- | --- |
| 1 | `--database` / `-d` |
| 2 | `GVCLASS_DB` environment variable |
| 3 | `database.path` in the configuration |
| 4 | `resources/` in the GVClass repository |

Relative database paths are resolved from the repository directory.

## Put the resource cache in a writable directory

Compact bundles store labels and reference proteins in Parquet files. GVClass
creates the label tables and marker FASTA files needed for a run in
`.gvclass_cache/` under the database directory.

For a read-only shared database, set a writable cache path in the same shell or
job script:

```bash
export GVCLASS_RESOURCE_CACHE="$HOME/.cache/gvclass"
mkdir -p "$GVCLASS_RESOURCE_CACHE"
```

This variable overrides `database.cache_path` in the configuration. A relative
`database.cache_path` is resolved inside the database directory.

## Check for database updates

```bash
pixi run setup-db
```

Setup and normal runs compare the installed version with the available download
source. GVClass checks the latest linked Zenodo record when available, or uses
the configured download source.

An interactive run asks before replacing an older complete database.
Non-interactive runs update it automatically. Check the database before starting
a batch job if downloads are restricted on compute nodes. Record the installed
`DB_VERSION` with your results.

Download URL, version, and checksum settings are listed in the
[configuration reference](../reference/configuration.md). Keep these settings
consistent with the selected archive. For `EUK-pEVE` reference labels, see
[Taxonomy and classification](../explanation/taxonomy.md#viral-panels-and-peve-references).
