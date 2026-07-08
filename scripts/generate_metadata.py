# init
# ----
import polars as pl
from Bio import SeqIO

from scripts.utils import handle_gzip

# logging
# -------
from log_config import configure_logger
import os
logger_name = os.path.splitext(os.path.basename(__file__))[0]
logger      = configure_logger(logger_name)

# snakefile parameters
# --------------------

GENOME_LIST             = snakemake.input['genome_list']
TAXONOMY_MAP            = snakemake.params['taxonomy_map']
METADATA_OUTPUT         = snakemake.output['metadata']

# functions
# ---------

@handle_gzip(logger)
def count_nucleotides(fasta_path):
    """
    Counts the number of nucleotides in a FASTA file.
    """
    return sum(
        len(record.seq) 
        for record 
        in SeqIO.parse(fasta_path, "fasta")
        )


def generate_metadata(genome_list, taxonomy_map, metadata_output):
    """
    Generates metadata for metapackage benchmarking.
    Uses columns: 
        accession,
        genome_size,
        taxonomy
    """

    logger.info("Reading genome list...")
    with open(genome_list, 'r') as f:
        genomes = [line.strip().split('\t')[1] for line in f.readlines()]

    logger.info(f"Found {len(genomes)} genomes.")
    for genome in genomes:
        logger.debug(f"Genome: {genome}")

    logger.info("Counting genome sizes...")
    genome_sizes = [
        count_nucleotides(genome) 
        for genome 
        in genomes
        ]

    logger.info("Creating metadata dataframe...")
    df = pl.DataFrame({
        'accession': [os.path.basename(g).split('.')[0] for g in genomes],
        'genome_size': genome_sizes
        })

    logger.info("Reading taxonomy map...")
    tax_df = (
        pl.read_csv(
            taxonomy_map, 
            separator='\t', 
            has_header=True)
        .select([
            pl.col('acc').alias('accession'), 
            pl.col('parsed_taxon_string').alias('taxonomy')]))

    print(tax_df)
    print(df)

    logger.info("Merging metadata with taxonomy...")
    merged_df = df.join(tax_df, on='accession', how='left', coalesce=True)
    logger.info("Writing metadata to output...")
    merged_df.write_csv(
        metadata_output, 
        separator='\t')

# main
# ----
def main():
    generate_metadata(GENOME_LIST, TAXONOMY_MAP, METADATA_OUTPUT)

# run
# ---
if __name__ == "__main__":
    main()