## Plan: Expand Benchmarking Workflow

Refactor `benchmarking_workflow` into a config-driven benchmark harness aligned with the current `metapackage_workflow` outputs. Preserve the existing simple simulation as a regression dataset, add independent complex simulations, scaffold accession-based real datasets with Kingfisher + Wagtail + taxonomy remapping, and introduce a tool adapter/evaluation layer with separate 18S and CoverM truth tracks.

**Current assessment**
- The current config enables only `coverage_definitions/simple`; the Snakefile derives `euk0...eukN` from files in one coverage directory.
- Only `euk-singlem` is executed. Existing converter scripts/environments show prior support for CoverM, Kraken/Bracken, Kaiju, MetaBuli, MetaPhlAn, mOTUs, Sourmash, and MAP2B, but they are not wired into this DAG.
- Simulation truth is generated from assigned ART genome coverages and the taxonomy map, then evaluated through Biobox/OPAL. There is no real-read acquisition, 18S truth, CoverM execution, or aggregate benchmark summary.
- `select_genomes` currently writes one FASTA path per line, while `generate_community.py` reads two tab-separated fields (`accession`, `fasta_path`): this is a blocking manifest contract bug.
- Coverage values are silently truncated to the number of genomes, and there is no validation that coverage rows, genome IDs, taxonomy, and selected genomes agree.
- All samples share one selected-genome manifest, so simple/complex datasets cannot currently have independent communities.
- Benchmark config points to the superseded `outputs/euk_packages/final_metapackage.smpkg` and `outputs/genomes/*/*.fasta.gz` conventions. Current authoritative paths are `metapackage_workflow/outputs/metadata/manifests/genomes_to_process.paths`, `metapackage_workflow/outputs/metadata/taxonomy/taxonomy_mapping.tsv`, `metapackage_workflow/outputs/genomes/{group|ncbi|p10k}/{genome_acc}.fasta.gz`, and `metapackage_workflow/outputs/singlem/metapackages/final.smpkg`.
- The current final-output path is nested as `outputs/{tool}/{tool}/{sample}...`; thread-specific results are deleted by the comparison rule, and there is no durable per-dataset/per-truth result index.

**Steps**

### Phase 1: Stabilise the existing benchmark
1. Make the benchmark consume the current metapackage workflow paths through one explicit path/config layer, preferably resolving paths relative to the workspace and allowing overrides for external assets. Reuse `metapackage_workflow/config/paths.yaml` conventions rather than duplicating legacy paths.
2. Define a canonical genome manifest schema with accession, FASTA path, source/group, and taxonomy. Fix selection and metadata generation to emit/read that schema consistently, derive accession from the manifest rather than filename heuristics, and fail on missing taxonomy or duplicate IDs.
3. Add input validation for coverage definitions: exact/declared genome count, non-negative finite coverages, required taxonomy, and explicit handling of zero-coverage genomes. Keep the current simple coverage files as regression fixtures.
4. Preserve thread benchmarking without deleting raw thread-specific outputs. Store final selected output separately and record whether profiles/OTU tables are identical; make the output layout `{dataset}/{tool}/{sample}/{threads}` plus a final result directory.
5. Add a dry-run/parse validation target and a small fixture test for the manifest, coverage, taxonomy, and path contracts before adding new datasets.

### Phase 2: Explicit simulated datasets
6. Replace the single coverage-folder/sample derivation with explicit dataset records, each declaring `id`, `kind: simulated`, genome-selection policy/seed, coverage source or generator, read simulator parameters, truth outputs, and reference/metapackage variant. Support separate manifests and seeds for every community; allow simple and complex datasets to share a policy only when explicitly configured.
7. Add `coverage_definitions/complex` as a first-class dataset family and add configurable realistic distribution generators (for example sparse log-normal/long-tail, uneven dominant taxa, and zero/low abundance cases). Materialise generated coverage tables so simulations are reproducible and inspectable.
8. Generalise the simulation rule to consume a dataset manifest and emit paired reads, genomewise truth, taxonomy-condensed truth, and an explicit truth metadata record. Keep ART parameters in config and avoid embedding `euk` naming assumptions in rule logic.

### Phase 3: Novelty-aware reference variants
9. Add a reference-variant abstraction for in-reference, partially novel, and out-of-reference communities. A dataset should declare which genomes are present in the simulated community and which are visible to each profiler/metapackage.
10. Implement a custom metapackage preparation step that filters/removes selected lineages and their package entries consistently at all taxonomy ranks, emits an exclusion/report manifest, and creates a dataset-specific metapackage. Do not mutate the current final metapackage in place. Define whether novelty means absent genomes, absent species/lineages, or both, and report the exact held-out IDs/taxa.
11. Add truth labels separating detected reference taxa, held-out taxa, and contamination/off-target taxa so novelty performance can be evaluated rather than collapsed into one OPAL score.

### Phase 4: Real characterized metagenome scaffold
12. Add a `real` dataset schema based on accession manifests rather than local FASTQ paths. Include accession, sample ID, read layout, download status, expected files, metadata, and optional reference database requirements.
13. Scaffold Kingfisher-download rules using the repository/project URL and pinned environment/version. The rules should be resumable and produce a manifest/checksum marker; keep accessions/configuration separate from implementation.
14. Scaffold Wagtail 18S processing as an independent branch: downloaded reads -> Wagtail outputs -> taxonomy remapping to the benchmark/metapackage taxonomy -> normalized 18S truth table. Keep the remapping script explicit about rank mapping, unresolved taxa, synonyms, and taxa absent from the metapackage.
15. Add an optional CoverM branch for real samples, with reference FASTA/index configuration and a normalized genome/taxon coverage truth table. Keep 18S-derived and CoverM-derived truth tracks separate; do not force them into one ground truth when they measure different biological signals.
16. Mark the real-data rules as scaffolded/unvalidated until accession examples, Wagtail command/version, and CoverM reference inputs are supplied. Provide config placeholders and clear failure messages rather than fake defaults.

### Phase 5: Tool comparison layer
17. Introduce a tool registry/config contract describing command, environment, database/reference inputs, read mode, output format, taxonomy mapping, and applicable dataset kinds. Keep SingleM as the first working adapter.
18. Prioritise adapters for Kraken (and Bracken if used), MetaPhlAn4, EukDetect, FunOMIC, and MiCoP. Add CoverM as a separate coverage estimator. Reuse existing converter scripts where their format/taxonomy assumptions match; repair or replace converters that assume bacterial/archaeal taxonomies, hard-coded ranks, or NCBI-only IDs.
19. Add eukaryote/multi-domain profiler candidates only behind explicit config entries, with a compatibility matrix recording whether each tool supports simulated genomes, real reads, novelty references, and 18S evaluation.
20. Make each tool emit a canonical condensed profile plus provenance (tool version, database, command, reference variant, input sample). Route all canonical profiles through common Biobox/OPAL conversion where valid, while retaining tool-native outputs.

### Phase 6: Evaluation and reporting
21. Parameterise OPAL/evaluation by dataset and truth track. For simulated data, evaluate against assigned genome/taxonomy truth; for real data, run separate 18S and CoverM comparisons and label them accordingly.
22. Add an aggregate report (TSV/JSON and optionally notebook figures) covering per-tool/per-dataset/per-thread runtime, memory, completeness/false positives, taxonomic rank metrics, novelty detection, and truth-track agreement. Do not compare 18S relative abundance directly with whole-metagenome coverage without a documented normalization/interpretation.
23. Update README/config documentation with the new path conventions, dataset schema, tool registry, accession workflow, novelty semantics, environment setup, and expected outputs.

**Relevant files**
- `benchmarking_workflow/snakefile` — replace implicit coverage-folder/sample logic with dataset-driven targets; fix manifest contract; add real-data, reference-variant, tool, truth-track, and reporting rules.
- `benchmarking_workflow/config.yaml` — current legacy inputs, thread/tool lists, resources, and path corrections; split into authoritative paths plus dataset/tool configuration as needed.
- `benchmarking_workflow/scripts/generate_community.py` — canonical manifest/coverage validation, independent dataset inputs, reproducible simulation, and explicit truth metadata.
- `benchmarking_workflow/scripts/generate_metadata.py` — consume the canonical manifest and current taxonomy mapping without filename-derived accession assumptions.
- `benchmarking_workflow/coverage_definitions/simple/*` and `coverage_definitions/complex/*` — regression and complex simulation definitions.
- `benchmarking_workflow/conversion_scripts/*` — adapter normalization; especially existing CoverM, Kraken, and MetaPhlAn converters and their taxonomy assumptions.
- `benchmarking_workflow/envs/*` — pin/add Kingfisher, Wagtail, Kraken, MetaPhlAn4, EukDetect, FunOMIC, MiCoP, CoverM, and any selected databases.
- `benchmarking_workflow/README.md` — document execution, dataset/tool configuration, current output layout, real-data scaffold, and limitations.
- `metapackage_workflow/config/paths.yaml` and `metapackage_workflow/OUTPUT_LAYOUT.md` — authoritative current path contracts to reuse/reference.
- `metapackage_workflow/outputs/metadata/manifests/genomes_to_process.paths`, `outputs/metadata/taxonomy/taxonomy_mapping.tsv`, and `outputs/singlem/metapackages/final.smpkg` — current benchmark inputs, subject to existence checks at runtime.

**Verification**
1. Run Snakemake parse/dry-run checks for the legacy simple dataset against the current metapackage paths; verify all expected inputs resolve and no old `outputs/euk_packages` path remains.
2. Run the manifest/coverage validation fixture and a one-sample simple simulation; verify paired FASTQ files, genomewise truth, condensed truth, taxonomy IDs, and reproducibility from a fixed seed.
3. Run one complex simulation with an independent genome manifest and confirm it does not reuse the simple dataset’s selected genomes or coverage file.
4. Run a novelty fixture with a dataset-specific filtered metapackage; verify held-out genomes/taxa are absent from the reference and retained in truth labels.
5. For each tool adapter, run a minimal fixture or command dry-run and verify canonical condensed output, provenance, taxonomy remapping, and failure behavior when its database is not configured.
6. Run the real-data scaffold in dry-run mode and verify Kingfisher, Wagtail, taxonomy-remap, and CoverM checkpoints form separate branches with explicit placeholder errors for missing accessions/inputs.
7. Validate OPAL outputs and aggregate reports across at least one simulated truth track and confirm separate 18S/CoverM real-data reports are not merged.

**Decisions**
- Real data will be scaffolded for future accession input rather than selecting public accessions now.
- 18S and CoverM are separate truth tracks; neither is silently treated as universally authoritative.
- Communities are independently generated with semi-random genome selection and reproducible seeds.
- Novelty support is included as a first-class reference-variant feature, with dataset-specific metapackages and explicit held-out truth.
- Initial comparison focus: Kraken, MetaPhlAn4, EukDetect, FunOMIC, MiCoP, plus CoverM; other eukaryote/multi-domain tools remain optional registry entries.
- No implementation edits are included in this planning pass.

**Further Considerations**
1. Before implementing real-data commands, pin the exact Kingfisher and Wagtail versions and provide one accession/read-layout example; the scaffold should otherwise remain config-driven.
2. Decide whether the initial benchmark metric should be OPAL-only or also include rank-aware precision/recall and abundance error; recommendation: add rank-aware metrics because novelty and separate truth tracks exceed OPAL’s single comparison model.
3. Confirm whether the current SingleM metapackage should be treated as the benchmark reference artifact or regenerated per novelty scenario; recommendation: preserve the current artifact and derive immutable scenario-specific copies.
