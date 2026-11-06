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
    manifest_df = pl.read_csv(
        genome_list,
        separator='\t',
        has_header=False,
        new_columns=['accession', 'fasta_path'],
    )
    genomes = manifest_df['fasta_path'].to_list()
    accessions = manifest_df['accession'].to_list()

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
        'accession': accessions,
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

    logger.info("Merging metadata with taxonomy...")
    merged_df = df.join(tax_df, on='accession', how='left', coalesce=True)
    missing_taxonomy = merged_df.filter(pl.col('taxonomy').is_null())
    if len(missing_taxonomy) > 0:
        missing = missing_taxonomy['accession'].to_list()
        raise ValueError(f"Missing taxonomy for selected genomes: {missing}")
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