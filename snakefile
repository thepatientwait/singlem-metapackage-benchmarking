
from pathlib import Path
# config

import sys
from pathlib import Path

CONFIG_PATH       = "config.yaml"
SCRIPTS_DIR       = "scripts"
METAPACKAGE_ROOT  = Path(workflow.basedir).parent / "metapackage_workflow"
PATHS_CONFIG_PATH = METAPACKAGE_ROOT / "config" / "paths.yaml"

from helpers import (
    create_resource_profile_factory,
)

def ensure_list(value):
    if isinstance(value, (list, tuple, set)):
        return list(value)
    return [value]

configfile: str(CONFIG_PATH)
configfile: str(PATHS_CONFIG_PATH)


# resources
BASE_MEM_MB    = config["resources"]["base_mem_mb"]
BASE_RUNTIME   = config["resources"]["base_runtime"]
TOTAL_ATTEMPTS = config["resources"]["total_attempts"]

RESOURCE_PROFILE   = create_resource_profile_factory(BASE_MEM_MB, BASE_RUNTIME, TOTAL_ATTEMPTS)
DEFAULT_RESOURCES  = RESOURCE_PROFILE()

THREAD_COUNTS               = [int(t) for t in ensure_list(config['benchmarking_threads'])]
TOOLS                       = config['tools']
# dirs
BENCHMARK_DIR               = Path(workflow.basedir) / config['benchmark_dir']
OUTPUT_DIR                  = Path(workflow.basedir) / config['output_dir']
METADATA_DIR                = Path(workflow.basedir) / config['metadata_dir']
FASTQ_DIR                   = OUTPUT_DIR / 'simulated_reads'
TRUTH_DIR                   = OUTPUT_DIR / 'truths'

# Ensure supporting folders exist (created at parse time if missing)
for _d in (METADATA_DIR, OUTPUT_DIR, BENCHMARK_DIR, FASTQ_DIR, TRUTH_DIR):
    _d.mkdir(parents=True, exist_ok=True)

# inputs from metapackage_workflow/config/paths.yaml
UPSTREAM_GENOME_MANIFEST    = METAPACKAGE_ROOT / config['download']['manifest']['genome']
SINGLEM_METAPACKAGE         = METAPACKAGE_ROOT / config['generate']['singlem']['metapackages']['final']
TAXONOMY_MAP                = METAPACKAGE_ROOT / config['download']['metadata']['taxonomy']['mapping']

# TODO: update to parse multiple coverage definiton classes
COVERAGE_DEFINITIONS_FOLDER = Path(workflow.basedir) / config['coverage_definitions_folder'][0]
GENOME_NUMBER               = max((sum(1 for _ in p.open('r')) for p in COVERAGE_DEFINITIONS_FOLDER.iterdir() if p.is_file()), default=0)
num_coverage_files          = len([p for p in COVERAGE_DEFINITIONS_FOLDER.iterdir() if p.is_file() and p.name.startswith('coverage')])
# num_coverage_files          = 1
DATASETS                    = list(['euk'+str(i) for i in range(num_coverage_files)])

# outputs
GENOME_FASTA_PATHS          = METADATA_DIR / config['genome_manifest']
BENCHMARK_INPUTS_READY      = METADATA_DIR / 'benchmark_inputs.ready'
OPAL_REPORTS                = OUTPUT_DIR / "{tool}" / "opal" / "{sample}.opal_report"


wildcard_constraints:
    threads="[0-9]+",
    tool="[a-zA-Z0-9_-]+",
    sample="euk[0-9]+"



rule all:
    input:
        str(BENCHMARK_INPUTS_READY),
        expand(str(OPAL_REPORTS), sample=DATASETS, tool=TOOLS)


rule check_benchmark_inputs:
    input:
        genome_manifest = str(UPSTREAM_GENOME_MANIFEST),
        taxonomy_map    = str(TAXONOMY_MAP),
    params:
        metapackage     = str(SINGLEM_METAPACKAGE),
    output:
        ready = str(BENCHMARK_INPUTS_READY),
    run:
        required_inputs = {
            'genome manifest': input.genome_manifest,
            'taxonomy map': input.taxonomy_map,
            'final SingleM metapackage': params.metapackage,
        }
        missing = [name for name, path in required_inputs.items() if not Path(path).exists()]
        if missing:
            raise FileNotFoundError(
                'Benchmark inputs are not ready; missing: ' + ', '.join(missing) +
                '. Complete metapackage_workflow first.'
            )
        Path(output.ready).write_text('Benchmark inputs validated\n', encoding='utf-8')


rule select_genomes:
    input:
        upstream_manifest = str(UPSTREAM_GENOME_MANIFEST),
    params:
        genome_count          = GENOME_NUMBER,
        random_seed            = 42,
    output:
        genome_list         = GENOME_FASTA_PATHS,
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    run:
        from numpy import random

        def accession_from_path(path):
            name = Path(path).name
            for suffix in ('.fasta.gz', '.fna.gz', '.fa.gz', '.fasta', '.fna', '.fa'):
                if name.endswith(suffix):
                    return name[:-len(suffix)]
            raise ValueError(f"Unsupported genome filename: {path}")

        random.seed(params.random_seed)
        with open(input.upstream_manifest, 'r', encoding='utf-8') as manifest:
            fasta_paths = []
            for line in manifest:
                raw_path = line.strip()
                if not raw_path:
                    continue
                fasta_path = Path(raw_path)
                if not fasta_path.is_absolute():
                    fasta_path = METAPACKAGE_ROOT / fasta_path
                fasta_path = fasta_path.resolve()
                if not fasta_path.is_file():
                    raise FileNotFoundError(
                        f"Genome listed by upstream manifest does not exist: {fasta_path}"
                    )
                fasta_paths.append(str(fasta_path))
        if len(fasta_paths) < params.genome_count:
            raise ValueError(
                f"Upstream manifest contains {len(fasta_paths)} genomes, "
                f"but {params.genome_count} are required"
            )
        selected_genomes = random.choice(fasta_paths, size=params.genome_count, replace=False)
        with open(output.genome_list, 'w', encoding='utf-8') as outfile:
            print(f"Selecting genomes into {output.genome_list}")
            for fasta_path in selected_genomes:
                outfile.write(f"{accession_from_path(fasta_path)}\t{fasta_path}\n")

rule generate_metadata:
    input:
        genome_list         = GENOME_FASTA_PATHS,
    output:
        metadata            = str(METADATA_DIR / 'genome_metadata.tsv'),
    params:
        taxonomy_map        = TAXONOMY_MAP,
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    script:
        'scripts/generate_metadata.py'


rule generate_community_and_reads:
    input:
        coverage_file       = lambda wildcards: str(COVERAGE_DEFINITIONS_FOLDER / f"coverage{wildcards.sample.replace('euk', '')}.tsv"),
        genome_list         = GENOME_FASTA_PATHS,
        metadata            = str(METADATA_DIR / 'genome_metadata.tsv'),
    output:
        read1               = str(FASTQ_DIR / "{sample}.1.fq.gz"),
        read2               = str(FASTQ_DIR / "{sample}.2.fq.gz"),
        condensed           = str(TRUTH_DIR / "{sample}.condensed"),
        genomewise_coverage = str(TRUTH_DIR / "{sample}.genomewise.csv"),
    params:
        art_binary          = 'art_illumina',
    threads: 10
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    conda:
        'envs/art.yml'
    script:
        'scripts/generate_community.py'

###############################################################################################
###############################################################################################

rule euk_singlem_run_to_profile:
    input:
        r1                  = str(FASTQ_DIR / "{sample}.1.fq.gz"),
        r2                  = str(FASTQ_DIR / "{sample}.2.fq.gz"),
        benchmark_inputs    = str(BENCHMARK_INPUTS_READY),
    params:
        db                  = str(SINGLEM_METAPACKAGE),
    benchmark:
        str(BENCHMARK_DIR / "{tool}/{sample}-{threads}threads.benchmark")
    output:
        report              = str(OUTPUT_DIR / "{tool}/{sample}-{threads}.profile"),
        otu_table           = str(OUTPUT_DIR / "{tool}/{sample}-{threads}.otu_table"),
    conda:
        "../singlem/extras/lyrebird_metapackage_creation/envs/singlem.yml"
    threads: 
        lambda wildcards: int(wildcards.threads)
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    log:
        str(OUTPUT_DIR / "logs/{tool}/{sample}-{threads}.log")
    shell:
        "singlem pipe "
        "--threads {threads} "
        "-1 {input.r1} "
        "-2 {input.r2} "
        "-p {output.report} "
        "--otu-table {output.otu_table} "
        "--output-extras "
        "--metapackage {params.db} &> {log}"


rule compare_and_cleanup_threading_results:
    input:
        profiles = lambda wildcards: [str(OUTPUT_DIR / f"{wildcards.tool}/{wildcards.sample}-{threads}.profile") for threads in THREAD_COUNTS],
        otu_tables = lambda wildcards: [str(OUTPUT_DIR / f"{wildcards.tool}/{wildcards.sample}-{threads}.otu_table") for threads in THREAD_COUNTS]
    output:
        final_profile = str(OUTPUT_DIR / "{tool}/{tool}/{sample}.profile"),
        final_otu_table = str(OUTPUT_DIR / "{tool}/{tool}/{sample}.otu_table")
    params:
        thread_counts = THREAD_COUNTS
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    shell:
        """
        # Compare all profile files
        FIRST_PROFILE="{input.profiles[0]}"
        FIRST_OTU="{input.otu_tables[0]}"
        ALL_IDENTICAL=true
        
        for PROFILE in {input.profiles}; do
            if ! diff "$FIRST_PROFILE" "$PROFILE" >/dev/null 2>&1; then
                echo "Profile files differ: $FIRST_PROFILE vs $PROFILE"
                ALL_IDENTICAL=false
                break
            fi
        done
        
        for OTU in {input.otu_tables}; do
            if ! diff "$FIRST_OTU" "$OTU" >/dev/null 2>&1; then
                echo "OTU table files differ: $FIRST_OTU vs $OTU"
                ALL_IDENTICAL=false
                break
            fi
        done
        
        if [ "$ALL_IDENTICAL" = true ]; then
            echo "All threading outputs are identical, using first one"
        else
            echo "Threading outputs differ, using first one anyway"
        fi
        
        # Copy first files to final output
        cp "$FIRST_PROFILE" {output.final_profile}
        cp "$FIRST_OTU" {output.final_otu_table}
        
        # Remove all threading-specific files
        rm {input.profiles} {input.otu_tables}
        """


rule truth_condensed_to_biobox:
    input:
        condensed           = str(TRUTH_DIR / "{sample}.condensed"),
    output:
        biobox              = str(TRUTH_DIR / "{sample}.condensed.biobox"),
    conda:
        "envs/singlem.yml"    
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    shell:
        "{workflow.basedir}/../bin/condensed_profile_to_biobox.py "
        "--input-condensed-table {input.condensed} "
        "--output-biobox {output.biobox}"

rule tool_condensed_to_biobox:
    input:
        profile             = str(OUTPUT_DIR / "{tool}/{tool}/{sample}.profile"),
        truth               = str(TRUTH_DIR / "{sample}.condensed.biobox"),
    output:
        biobox              = str(OUTPUT_DIR / "{tool}/biobox/{sample}.biobox"),
    conda:
        "envs/singlem.yml"
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    shell:
        "{workflow.basedir}/../bin/condensed_profile_to_biobox.py "
        "--input-condensed-table {input.profile} "
        "--output-biobox {output.biobox} "
        "--template-biobox {input.truth} "

rule opal:
    input:
        biobox              = str(OUTPUT_DIR / "{tool}/biobox/{sample}.biobox"),
        truth               = str(TRUTH_DIR / "{sample}.condensed.biobox"),
    params:
        output_dir          = str(OUTPUT_DIR / "{tool}"),
        output_opal_dir     = str(OUTPUT_DIR / "{tool}/opal/{sample}.opal_output_directory"),
    output:
        report              = str(OUTPUT_DIR / "{tool}/opal/{sample}.opal_report"),
    conda:
        'envs/opal.yml'
    retries: TOTAL_ATTEMPTS
    resources: **DEFAULT_RESOURCES
    shell:
        "opal.py -g {input.truth} -o {params.output_opal_dir} {input.biobox} || echo 'expected opal non-zero existatus'; mv {params.output_opal_dir}/results.tsv {output.report} && rm -rf {params.output_opal_dir}"

