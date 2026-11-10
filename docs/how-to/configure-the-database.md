# Configure the database

Run these commands from the GVClass repository to install or select a reference
database. Database and software versions are separate.

## Install in the default location

```bash
pixi run setup-db
cat resources/DB_VERSION
```

Setup downloads and verifies the archive in a private staging directory next to
the database, extracts it there and checks that every required file is present.
Only then does it swap the new tree into place. A failed download, a checksum
mismatch or a missing file leaves an installed database untouched, and setup
removes the staging directory in every case. Setup never reads or writes a file
named like the archive next to the database. The bundle is about 1.5 GB and is
published on [Zenodo](https://doi.org/10.5281/zenodo.23381313) as 24 parts of 64 MiB.
Setup checks each part and the joined archive. A network error, an HTTP 5xx
answer or a checksum mismatch repeats the download of that part up to three times,
after 2, 8 and 30 seconds. Each request times out after 60 seconds. `DB_VERSION`
records the installed version.
Setup preserves this file exactly as archived and checks that its version matches
the selected source before replacing an existing database. An archive without
`DB_VERSION` receives a marker containing the source version and a final newline.

Setup holds an exclusive lock on `.<database name>.lock` in the parent directory
of the database while it installs or updates. Array jobs that start together are
safe: the second process waits, checks the database again and skips the download
when the first process has installed it.

Setup uses one exact source. An explicit environment URL has precedence over the
configuration block, which has precedence over the built-in v2.1.0 source. GVClass
does not follow Zenodo's latest release or fall back to an older bundle. See the
[configuration reference](../reference/configuration.md#database-keys) for the
environment variables and checksum rules.

## Install the contig-screen add-on

The optional contig-screen add-on installs next to the database. Install it with the
database, or afterwards. A database update keeps an installed `contig_screen/`
directory. The add-on installer takes its own lock, `<database>/.contig_screen.lock`,
inside the database directory, so it works when only the database directory is
writable. A database install or update takes that lock before it moves an installed add-on.

```bash
pixi run setup-db --with-screen
./gvclass screen --install
```

GVClass reads the add-on from `<database>/contig_screen/`. In the Apptainer image,
which has a read-only database, the add-on sits in a host directory that the image
mounts into the database. See
[Screen bins for contaminating contigs](screen-contigs.md).

## Use a shared database

1. Set `GVCLASS_DB` to an absolute path. Replace the example path with a directory
   on your shared filesystem.

    ```bash
    export GVCLASS_DB=/path/to/shared/gvclass_resources
    ```

2. Set up the database at that location. If a complete copy exists, setup
   reuses it.

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
`.gvclass_cache/` under the database directory, in one subdirectory per resource
identity. A new resource version starts a new subdirectory; old ones stay until
you remove them. Validation accepts reference proteins from `parquet/faa.parquet`
or from nonempty legacy `database/faa/*.faa` files.

For a read-only shared database, set a writable cache path in the same shell or
job script:

```bash
export GVCLASS_RESOURCE_CACHE="$HOME/.cache/gvclass"
mkdir -p "$GVCLASS_RESOURCE_CACHE"
```

This variable overrides `database.cache_path` in the configuration. A relative
`database.cache_path` is resolved inside the database directory.

## Check for database updates

GVClass never replaces a complete database without your consent. Normal runs and
`pixi run setup-db` compare the installed `DB_VERSION` with the effective download
source: the environment override, the configured source, or the built-in v2.1.0
default. A newer Zenodo release does not matter.

| Situation | Behavior |
| --- | --- |
| Database missing or incomplete | Setup installs the pinned source. |
| Complete, installed version equals or exceeds the pinned one | Nothing happens. |
| Complete, older, in a terminal | GVClass asks whether to update. The default is No. |
| Complete, older, not in a terminal | GVClass prints one warning with both versions and the update command, then runs with the installed database. |
| Complete, `DB_VERSION` missing or not a version number | GVClass treats it as a custom database, prints one line and never updates it. |

To update on purpose, run:

```bash
pixi run setup-db --update
```

The update replaces a complete database with the pinned source, whatever the
installed version. Check the database before you start a batch job if downloads are
restricted on compute nodes. Record the installed `DB_VERSION` with your results.

Download URL, version, and checksum settings are listed in the
[configuration reference](../reference/configuration.md). Keep these settings
consistent with the selected archive. For `EUK-pEVE` reference labels, see
[Taxonomy and classification](../explanation/taxonomy.md#viral-panels-and-peve-references).
