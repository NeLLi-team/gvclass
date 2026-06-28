# Getting started

By the end of this page you will have installed GVClass with [Pixi](https://pixi.sh), downloaded its reference database once, classified three bundled example genomes, and read the summary table that GVClass writes. Every command below is copy-paste ready and runs the same example, so you can check your output against the results shown here.

!!! note "What you need"
    A Linux machine (`linux-64`), plus `git` and `curl`. Pixi installs Python and every other dependency into a project-local environment, so nothing else has to be set up first.

## 1. Install Pixi

Pixi manages the GVClass environment. Install it with the official script.

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

Open a new shell so the `pixi` command is on your `PATH`, then confirm it is available.

```bash
pixi --version
```

You should see a version string such as `pixi 0.40.0`. Once the command is found, Pixi is ready.

## 2. Clone GVClass and install dependencies

Clone the development branch and move into the repository.

```bash
git clone -b gvclass-dev https://github.com/NeLLi-team/gvclass.git
cd gvclass
```

Install the environment. Pixi reads `pixi.toml`, resolves every pinned tool (pyrodigal, pyhmmer, VeryFastTree, and the rest), and builds a local environment.

```bash
pixi install
```

The first install downloads packages and takes a few minutes. When it finishes, the environment is locked and reused on every later command. Run GVClass commands from this `gvclass` directory so the launcher can find `src/`.

## 3. Download the database

GVClass classifies against a reference database of giant virus orthologous groups. Download it once.

```bash
pixi run setup-db
```

!!! note "One download per database location"
    This fetches the v1.7.1 resource bundle (about 2 GB) and verifies its checksum before unpacking. It runs once per database location; later runs reuse the files in `resources/`. To point GVClass at a shared copy, see [Configure the database](../how-to/configure-the-database.md).

When the download and checksum verification finish, a `resources/` directory holds the database.

## 4. Classify the example genomes

GVClass ships with three inputs in `example/`: two nucleotide bins (`.fna`) and one protein set (`.faa`). Run the bundled example.

```bash
pixi run example
```

This classifies the three files in `example/` into `example_results/` using 8 threads, with the default VeryFastTree tree method, fast mode, and sensitive mode all on. To change any of these, see [Tune speed and accuracy](../how-to/tune-speed-and-accuracy.md) and the [CLI reference](../reference/cli.md).

The run prints a configuration banner, per-query progress, and a completion message. Your output should look something like the block below (the layout is representative; the exact numbers will differ).

```text
    ================================================================
    GVClass - Giant Virus Classification Tool
    ================================================================

============================================================
                    Pipeline Configuration
============================================================
Config file: config/gvclass_config.yaml
Query directory: /home/you/gvclass/example
Output directory: /home/you/gvclass/example_results
Database: /home/you/gvclass/resources
Threads: 6 used / 8 requested (Workers: 3 × 2 threads)
Tree method: veryfasttree
Fast mode: True
Sensitive mode: True
Resume mode: DISABLED
Total queries: 3
============================================================
Completeness mode: novelty-aware

Progress: [ 33%] Query 1/3 | Memory: 540MB (peak: 612MB)
Progress: [ 67%] Query 2/3 | Memory: 690MB (peak: 720MB)
Progress: [100%] Query 3/3 | Memory: 705MB (peak: 800MB)

Pipeline completed successfully!
Validating pipeline outputs...
All outputs validated successfully
Generating combined summary...
Combined summary written to: example_results/gvclass_summary.tsv
CSV summary written to: example_results/gvclass_summary.csv
```

Notice that an `example_results/` directory now exists. The combined table you will read next is `example_results/gvclass_summary.tsv`.

## 5. Read your first result

Open `example_results/gvclass_summary.tsv`. It has one row per query and 44 columns. A few of those columns tell you most of what you want at a glance.

| query | taxonomy_majority | taxonomy_confidence | estimated_completeness | estimated_contamination | contamination_type | ttable |
|-------|-------------------|---------------------|------------------------|-------------------------|--------------------|--------|
| AC3300027503___Ga0255182_1000024 | d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;o_Imitervirales;f_IM_01;g_g2284;s_S392 | high | 100.00 | 0.00 | clean | codemeta |
| GVMAG-S-1096109-37 | d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;o_Pimascovirales;f_PM_01;g_g94;s_S679 | high | 82.86 | 0.08 | clean | codemeta |
| PkV-RF01 | d_NCLDV;p_Nucleocytoviricota;c_Megaviricetes;o_Imitervirales;f_IM_19;g_g1787;s_singleton | high | 100.00 | 0.00 | clean | no_fna |

`taxonomy_majority` is the full lineage, from domain (`d_`) through species (`s_`), produced by the per-marker single-gene-tree nearest-neighbor majority vote. All three queries classify with `high` confidence and a `clean` contamination type, which is what a high-quality giant virus genome looks like. For what completeness and contamination mean, see [Completeness and contamination](../explanation/quality-metrics.md).

The summary also reports one column per rank (`domain`, `phylum`, `class`, `order`, `family`, `genus`, `species`), but these do not simply repeat the lineage. Each holds the vote behind that rank as a distribution of single-protein-tree calls with counts and percentages. For `AC3300027503___Ga0255182_1000024` the `order` column reads `NCLDV__Imitervirales:19(95.00%),PLV_other:1(5.00%)`: nineteen of the twenty single-protein trees placed the query in Imitervirales and one landed elsewhere, so Imitervirales becomes the order in `taxonomy_majority`. See [Output files and columns](../reference/output.md) for every column.

`PkV-RF01` was a protein (`.faa`) input, so GVClass skipped gene calling: its genetic code shows `no_fna`, and the GC content and coding-density columns are not computed for it. The two `.fna` bins were gene-called, so they report a real genetic code (`codemeta`) and nucleotide statistics.

That is a full GVClass run. `example_results/` now holds the per-query files alongside the combined summary; exact vote counts and distances can shift slightly from run to run because tree inference is multithreaded.

## Next steps

- Classify your own data with [Classify a directory of bins](../how-to/classify-bins.md).
- Learn what every column means in [Output files and columns](../reference/output.md).
- Understand the method behind the numbers in [How GVClass works](../explanation/how-it-works.md).
