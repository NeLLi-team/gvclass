# Porting notes: cellular contig screen

The package is a port of the research prototype in
`tasks/contig-contamination-model-20261015/` (git-ignored). Each prototype module is a
module here with the same function names where the function was kept, so that a diff
against the prototype stays readable. Only the code needed to score a bin is ported:
training, cross-validation, hold-out and evaluation code is left out.

## Module map

| Prototype | Package | Kept | Left out |
| --- | --- | --- | --- |
| `predict_bin.py` | `cellular.py` | stage 1 stack, stage 2, DIAMOND searches, set logits, summary | the non-stacked path, `--set-dir`, `--set-variant`, `--context-dir`, `--threshold2-low` |
| `extract_genes.py` | `extract_genes.py` | `tile_contigs`, `call_genes`, `hmm_search`, `reference_search`, `write_tables` | `draw_windows`, `draw_contigs`, `Window`, the CLI, the done markers |
| `build_features.py` | `build_features.py` | `s_`, `c_`, `h_` features at hold-out `none` | all hold-out levels, `source_labels`, `fine` mode, the process pool |
| `ref_hits.py` | `ref_hits.py` | best hit per class, `piece_features`, best-hit table | all explicit and inferred hold-out code, caps, the lost-protein table |
| `extra_features.py` | `extra_features.py` | `Embedder`, `contig_proteins`, the `e_` features | k-mer frequencies (`k4_`, never used by the model), `draw_contigs` mode |
| `gene_embeddings.py` | `gene_embeddings.py` | embedding of all proteins, token extension | shard files, `--max-kb` subsets, `build`/`embed` CLI |
| `hit_tokens.py` | `hit_tokens.py` | token matrix, ablations | metrics (AUROC, ROC points), the CLI |
| `set_model.py` | `set_model.py` | `SetNet`, `Scaler`, `EmbeddingScale`, scoring | training, `cv`, `fit`, `nested`, folds |
| `stack_context.py` | `stack_context.py` | `stack_scores` | out-of-fold tables, nested tables, threshold fitting |
| `context_stage.py` | `context_stage.py` | `bin_context`, `context_table`, `tilings` | synthetic bins, training, `bin_thresholds` |
| `thresholds.py` | `thresholds.py` | `window_cover`, `sliding_calls` | `threshold_for_fpr`, hysteresis, length thresholds |
| `export_fasta.py`, `prepare_queries.py` | `cellular.py` (`query_index`, `write_queries`) | the protein FASTA and the per-window gene count | window FASTA files, gzip, the CLI |
| (new) | `layout.py`, `torch_stack.py` | manifest, hash checks, optional-dependency error | |

## Logic changes

Every change below keeps the numbers of the prototype. The equivalence run (see the
report of the port) compares the stage 1 and stage 2 scores of four bins.

| Module, line | Change | Reason |
| --- | --- | --- |
| `cellular.py:142` `main` | Prints an error and returns 2 for a missing torch or a manifest failure. | Design decision 3 and 4: clear refusal. |
| `cellular.py:202` `window_features` | One process, no subprocess; intermediate tables stay in memory except the three extraction tables. The window FASTA and `export_fasta.py` are gone; the proteins come from the gene calls of `extract_genes`. | Design decision 1 and 5. The prototype called the genes a second time with the same function and seeds; the tables are identical. |
| `cellular.py:202` `window_features` | DIAMOND and the set model always run. | The prototype ran them for `r_` features and `--set-dir`. The stack needs both. |
| `cellular.py:311` `diamond_search` | Binary found beside the interpreter, on `PATH`, or by `--diamond`; `--tmpdir` is the work directory; the search runs with captured output. | The prototype had a fixed path. |
| `cellular.py:404` `open_text` | `fasta_bp` and `contig_lengths` read gzipped files. | `read_fasta` already did; the summary failed on `.gz` input. |
| `extract_genes.py:54` `extract_bin` | Raises `ValueError` for a bin without a contig of 3 kb. | The prototype wrote a done marker and failed later at `pd.concat([])`. |
| `extract_genes.py:217` `load_models`, `:244` `reference_search` | The HMM catalog and the marker references are read from the verified bundle (`layout.py`), not from fixed paths. | Design decision 2 and 3. |
| `build_features.py:80` `load_labels`, `:93` `unit_features`, `:109` `annotate_hits` | Only hold-out `none`. The `r2` to `r7` lineage keys and the `proxy_*` columns are not built; `annotate_hits` takes the label table as an argument instead of a module global; the output keeps `holdout = "none"`. | The other levels and the keys feed only the training hold-outs. No model column uses them. |
| `ref_hits.py:56` `read_hits`, `:93` `collect_best` | An empty DIAMOND file gives no chunks and an empty best-hit table. | The prototype failed on an empty file. |
| `ref_hits.py:81` `best_per_class` | No exclusion, cap, or tally arguments. | With `--no-holdout` and no caps the prototype dropped nothing but references without a class. |
| `ref_hits.py:199` `reference_features` | The best-hit table (`best_table`, the former `write_best`) is returned, not written. | In memory instead of `best.parquet`. |
| `extra_features.py:124` `contig_proteins`, `:143` `window_embeddings` | The sample is drawn from the gene rows of `extract_genes` and not from a second pyrodigal call. Same filters, same hashes, same `rng.sample`. | One gene call per window. pyrodigal 3.3.0 and 3.7.1 gave the same genes on the three test bins. |
| `extra_features.py:95` `content_seed`, `:101` `draw_plan`; `cellular.py:202` `window_features` | The sampling seed is a digest of the contig sequences and the windows are drawn in content order, not file order. `--sample-seed NAME` restores the prototype seed and order. | The release score must not depend on the file name or contig order (lead's decision). |
| `extra_features.py:27` `Embedder` | Loads from a local directory with `local_files_only=True`. | No network at run time (design decision 2). |
| `gene_embeddings.py:49` `embed_proteins`, `:85` `extend_tokens` | Proteins come as `(name, sequence)` pairs; no shard files; the embeddings keep the float16 rounding of the stored shards. | In memory. |
| `hit_tokens.py:140` `window_tokens` | Replaces `main`: windows sorted by `contig_id`, label `y = 1`, same error for hits that match no gene. | In memory. |
| `set_model.py:176` `load_fitted` | Reads each file once, checks it against the manifest, then `torch.load(BytesIO, weights_only=True)`. | Design decision 3. |
| `set_model.py:200` `score_contigs` | Scores the requested variant only (`hits_struct`). | The prototype scored all saved variants and used one. |
| `stack_context.py:15` `stack_scores` | The clip value is the constant `PROBABILITY_CLIP`. | Same value, 1e-6. |
| `thresholds.py:28` `sliding_calls` | No `low` (hysteresis) argument. | `--threshold2-low` is not offered. |

## Seeded protein sample

The boosting model reads the mean ESM-2 embedding of at most 8 proteins per window.
The prototype seeded the draw with the file stem and drew in file order, so a renamed
or reordered copy of a bin scored differently (the prototype's runs A and B, the same
bin under two names, differ by up to 0.294 in `p_gbm` and 0.435 in `p_cellular`).

The port draws with `random.Random(f"{SEED}|esm|{digest}")`, where `digest` is the
SHA-256 of the sorted list of the SHA-256 hashes of the sequences (upper case) of the
contigs of 3 kb or more (`extra_features.content_seed`). The windows are drawn in the
order of (sequence hash of their contig, start) (`extra_features.draw_plan`). The
proteins of a window stay sorted by hash before `rng.sample`. The score therefore does
not depend on the file name, the contig names or the contig order. Two contigs with
the same sequence are the one case that still depends on file order.
`--sample-seed NAME` restores the prototype seed (`NAME` is the bin file stem) and the
file order of the windows; it exists to reproduce the numbers of the research
prototype. The set model embeds every protein and does not sample.

## Resources and integrity

`layout.py` documents the tree and its constants. Every file the screen reads is listed
in `contig_screen/manifest.json` (relative path, size, SHA-256; paths are relative to
the bundle root, so the bundle's own `hmm/combined.hmm` and
`parquet/labels/labels.parquet` are listed too). `MANIFEST_SHA256` in `layout.py` pins
the manifest. It is `None` until the release build sets it; with `None` the screen
refuses to run unless `--manifest-sha256` gives a digest. Joblib, torch, `.npz` and JSON
files are read once, hashed, and deserialized from those bytes. Files of 256 MiB or
more (the DIAMOND databases) are checked by size, or by hash with `--verify-all`.

The HMM catalog and the label table are read from the GVClass bundle root, not copied
into `contig_screen/`, so the release archive holds them once.

## Environment of the equivalence run

pixi environment `default` of this worktree (python 3.11): pandas 2.3.3, numpy 2.4.2,
scikit-learn 1.7.2, joblib 1.5.3, pyhmmer 0.12.0, pyrodigal 3.3.0, pyswrd 0.3.1,
biopython 1.85, pyarrow 22.0.0, scipy 1.17.0. The packages that the default
environment lacks came from `work/plm_env` of the prototype, appended to `sys.path`:
torch 2.14.1+cpu, transformers 5.19.0, tokenizers 0.23.3, safetensors 0.8.0,
huggingface_hub 2.2.0. DIAMOND 2.1.24 (the prototype's binary; the default environment
has 2.1.23). Release 2.1.0 adds `pytorch-cpu` and `transformers` to the default pixi
environment (`pixi.toml`, `pixi.lock`); the equivalence runs of the release use that
environment (torch 2.14.1, transformers 5.17.0), where the set logits differ from the
prototype's by up to 2e-6.

# Porting notes: giant-virus contig screen (model v9)

The detector of contigs of another giant virus in a giant-virus bin is a port of
`predict_gv_contamination.py` and its import closure in the same prototype folder. It
scores only. Training, cross-validation, synthetic bins, hold-outs and priors fitting are
left out.

## Module map

| Prototype | Package | Kept | Left out |
| --- | --- | --- | --- |
| `predict_gv_contamination.py` | `giant_virus.py` | `score_bin`, `report_union`, `main`, the printed summary | `--ref-dir`, `--model-dir`, `--cellular-low` |
| `gv_genes.py` | `gv_genes.py` | `read_fasta`, gene calls and catalog assignment (as `call_genes`) | the table CLI, the parquet files, the contig table (`gv_contigs`) |
| `gv_order.py` | `gv_order.py` | composition, codon and structure vectors, `Taxonomy`, `Contig`, `Priors`, `IdentityCache`, `bin_features`, `taxonomy_features` | `make_priors`, `priors_from_summaries`, `read_fasta` |
| `gv_order_build.py` | `gv_order_build.py` | `Genome` (codon counts, `contigs`), `add_taxa` | synthetic bins, folds, donors, priors fitting, the CLI |
| `gv_taxa.py` | `gv_taxa.py` | `best_catalog_hits`, `gene_table`, `taxa_for_bin` | fold and unit exclusions, the hold-out inputs, the CLI |
| `gv_proteins.py` | `gv_taxa.py` (`fasta_text`) | the FASTA text of the proteins | the second gene call and its comparison, the CLI |
| `ncldv_taxonomy_panel.py` | `gv_taxa.py` (`codes`) | family and genus code of a lineage | the panel evaluation |
| `thresholds.py` | `thresholds.py` (cellular port) | `sliding_calls` | hysteresis (`--cellular-low`) |
| (new) | `gv_model.py` | bundle and threshold loading | |

## Logic changes

None of the changes alters a score; the equivalence run compares `p_foreign` and the
calls of four bins with the prototype.

| Module, line | Change | Reason |
| --- | --- | --- |
| `gv_genes.py:65` `call_genes` | Replaces `process`. The gene table stays in memory; the HMM catalog comes from the verified bundle; the proteins of all genes are returned beside the table. | No parquet files, no done markers. The proteins serve the DIAMOND search. |
| `gv_taxa.py:76` `taxa_for_bin` | The search queries are the proteins of the same gene calls. The prototype called the genes a second time (`gv_proteins.called_proteins`) and checked that the calls matched row by row. | One gene call; the check is true by construction. |
| `gv_taxa.py:76` `taxa_for_bin` | DIAMOND runs through `cellular.diamond_binary` and `cellular.diamond_search` (same options, plus `--tmpdir`); the database and the catalog table are read from the verified bundle (`ref/ref7x.dmnd`, `ref/ref_genomes_x.tsv`). | Design decision 2 and 3. The prototype used a fixed binary and `--ref-dir`. |
| `gv_taxa.py:40` `best_catalog_hits` | No exclusion argument. | The prototype passed `{}` and no caps; `ref_hits.best_per_class` of the cellular port is that case. |
| `gv_taxa.py:52` `gene_table` | The output has no `genome_id` column; `row` is the key. | The table belongs to one bin and stays in memory. |
| `gv_order_build.py:28` `add_taxa` | Takes the gene table of `gv_taxa` as a frame, not the path of a parquet file with a per-genome cache. | In memory. The reindex and the fill values are those of the prototype. |
| `gv_order_build.py:45` `Genome` | Built from the contig sequences and the gene table; no file reads. | In memory. |
| `gv_order.py` `Contig`, `bin_features` | The `foreign` field and the `foreign` column are gone. | The training label; always 0 for a bin and not a model column. |
| `gv_model.py:35` `BundleUnpickler` | `model.joblib` is read once, hashed, and unpickled from those bytes; the prototype's class name `gv_order.Priors` maps to `gv_order.Priors` of this package. | The file holds a pickled `Priors` object. No module named `gv_order` has to exist, and no other project class is accepted. |
| `gv_model.py:43` `load_bundle` | A model without priors raises `ValueError`, which `main` reports with exit code 2. | The prototype exited through `SystemExit`. |
| `giant_virus.py:115` `main` | `--resources DIR` replaces `model dir` and `--ref-dir`. The threshold is key `0.01` of `giant_virus/thresholds.json` unless `--threshold` is given. Errors print `error: ...` to stderr and return 2. | Design decisions 1, 2 and 4. |
| `giant_virus.py:115` `main` | `--cellular-threshold` is optional: the default is the stage 2 threshold of the cellular `thresholds.json`, when the bundle lists it. Without either, the union is not computed (as in the prototype without `--cellular-threshold`). `--cellular-low` is not offered. | The cellular port has no hysteresis (`thresholds.sliding_calls`). |
| `giant_virus.py:115` `main` | The output directory is created. | |

## Resources

`layout.py` names the files: `contig_screen/giant_virus/model.joblib` and
`contig_screen/giant_virus/thresholds.json` (`GV_FILES` lists them with the three shared
files). The catalog DIAMOND database and table are the cellular screen's
`ref/ref7x.dmnd` and `ref/ref_genomes_x.tsv`, and the marker HMMs are the bundle's
`hmm/combined.hmm`; nothing is copied. The label table is not read. The model file holds
the fitted `HistGradientBoostingClassifier`, its 51 feature columns and the priors (how
often each marker is duplicated in single reference genomes, and the mean and spread of
the composition, codon and structure vectors); the prototype reads no other table.
`scripts/build_contig_screen_resources.py --gv-dir DIR` places the two files (`DIR`
holds `model.joblib` and `thresholds.json`) and lists them in the manifest. A bundle
that lists `giant_virus/model.joblib` must list the other four files as well.

## Environment of the equivalence run

pixi environment `default` of this worktree, the prototype's DIAMOND 2.1.24, the
prototype's models (`work/gv_order/v9`). The screen needs no torch.
