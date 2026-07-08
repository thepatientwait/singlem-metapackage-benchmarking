# init
# ----
import sys
import subprocess
import tempfile
import gzip
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import polars as pl

# Add parent directory to path for utils import
sys.path = [str(Path(__file__).parent.parent)] + sys.path

from scripts.utils import handle_gzip

# logging
# -------
from log_config import configure_logger
logger_name = Path(__file__).stem
logger = configure_logger(logger_name)

import logging
logger.setLevel(logging.DEBUG)

# snakefile parameters
# --------------------
COVERAGE_FILE               = snakemake.input['coverage_file']
GENOME_LIST                 = snakemake.input['genome_list']
METADATA_FILE               = snakemake.input['metadata']

ART_BINARY                  = snakemake.params['art_binary']

THREADS                     = snakemake.threads

OUTPUT_READ1                = snakemake.output['read1']
OUTPUT_READ2                = snakemake.output['read2']
OUTPUT_CONDENSED            = snakemake.output['condensed']
OUTPUT_GENOMEWISE_COVERAGE  = snakemake.output['genomewise_coverage']

# constants
# ---------
READ_LENGTH = 150
FRAGMENT_LENGTH = 400
FRAGMENT_STDEV = 10

# functions
# ---------

@handle_gzip(logger)
def load_coverage_data(coverage_file_path):
    """
    Load coverage values from file.
    
    Args:
        coverage_file_path: Path to coverage file with columns ['otu', 'coverage']
        
    Returns:
        List of coverage values
    """
    logger.info(f"Loading coverage data from {coverage_file_path}")
    
    try:
        coverage_df = pl.read_csv(
            coverage_file_path, 
            separator='\t', 
            has_header=False, 
            new_columns=['otu', 'coverage']
        )
        coverage_values = coverage_df['coverage'].to_list()
        
        logger.info(f"Loaded {len(coverage_values)} coverage values")
        logger.debug(f"Coverage values: {coverage_values}")
        
        return coverage_values
        
    except Exception as e:
        logger.error(f"Failed to load coverage data: {e}")
        raise


@handle_gzip(logger)
def load_genome_metadata(genome_list_path, metadata_file_path):
    """
    Load and merge genome list with metadata.
    
    Args:
        genome_list_path: Path to genome list file
        metadata_file_path: Path to metadata file with genome sizes and taxonomy
        
    Returns:
        polars.DataFrame: Merged genome metadata
    """
    logger.info("Loading genome list and metadata...")
    
    # Load genome list
    genome_df = pl.read_csv(
        genome_list_path, 
        separator='\t', 
        has_header=False, 
        new_columns=['accession', 'fasta_path']
    )
    logger.info(f"Loaded {len(genome_df)} genome entries")
    
    # Load metadata
    metadata_df = pl.read_csv(
        metadata_file_path, 
        separator='\t', 
        infer_schema_length=100000, 
        ignore_errors=True
    ).select('accession', 'genome_size', 'taxonomy')
    logger.info(f"Loaded {len(metadata_df)} metadata entries")
    
    # Merge datasets
    merged_df = genome_df.join(metadata_df, on='accession', how='inner')
    logger.info(f"Successfully merged {len(merged_df)} genome records")
    
    # Make fasta paths absolute using pathlib
    merged_df = merged_df.with_columns(
        pl.col('fasta_path').map_elements(lambda x: str(Path(x).resolve()), return_dtype=pl.Utf8)
    )
    
    return merged_df


def generate_simulation_commands(genome_metadata_df, coverage_values, art_binary_path, read_length, fragment_length, fragment_stdev, temp_dir):
    """
    Generate ART simulation commands for each genome.
    
    Args:
        genome_metadata_df: polars DataFrame with genome information
        coverage_values: List of coverage values to assign
        art_binary_path: Path to ART binary
        read_length: Length of simulated reads
        fragment_length: Mean fragment length
        fragment_stdev: Fragment length standard deviation
        temp_dir: Path object for temporary directory
        
    Returns:
        tuple: (simulation_commands, updated_genome_df)
    """
    logger.info("Generating simulation commands...")
    
    # Add coverage values to dataframe
    num_genomes = len(genome_metadata_df)
    coverage_subset = coverage_values[:num_genomes]
    
    genome_metadata_with_coverage = genome_metadata_df.with_columns(
        pl.Series('coverage', coverage_subset)
    )
    
    simulation_commands = []
    simulated_reads_dir = temp_dir / 'simulated_reads'
    
    for idx, row in enumerate(genome_metadata_with_coverage.rows(named=True)):
        fasta_path = row['fasta_path']
        coverage = row['coverage']
        
        # Use absolute path for output in temp directory
        output_prefix = simulated_reads_dir / f"{idx}-"
        
        command = (
            f"{art_binary_path} -ss HSXt -i {fasta_path} -p "
            f"-l {read_length} -f {coverage} -m {fragment_length} "
            f"-s {fragment_stdev} -o {output_prefix}"
        )
        logger.debug(f"Simulation command for genome {row['accession']}: {command}")
        simulation_commands.append(command)
    
    logger.info(f"Generated {len(simulation_commands)} simulation commands")
    return simulation_commands, genome_metadata_with_coverage


def write_output_files(genome_metadata_df, coverage_file_basename, condensed_output_path, genomewise_output_path):
    """
    Write condensed and genome-wise coverage output files.
    
    Args:
        genome_metadata_df: polars DataFrame with genome and coverage information
        coverage_file_basename: Basename of coverage file for sample naming
        condensed_output_path: Path for condensed output file
        genomewise_output_path: Path for genome-wise coverage file
    """
    logger.info("Writing output files...")
    
    # Write genome-wise coverage file
    genomewise_df = genome_metadata_df.select([
        pl.col('accession'),
        pl.col('coverage'), 
        pl.col('fasta_path').alias('fasta'),
        pl.col('taxonomy')
    ])
    
    genomewise_df.write_csv(
        genomewise_output_path,
        separator='\t'
    )
    
    # Calculate taxonomy-wise coverage totals using polars
    taxonomy_coverage_df = genome_metadata_df.group_by('taxonomy').agg(
        pl.col('coverage').sum().alias('total_coverage')
    ).with_columns(
        pl.lit(coverage_file_basename).alias('sample')
    ).select([
        pl.col('sample'),
        pl.col('total_coverage').alias('coverage'),
        pl.col('taxonomy')
    ])
    
    # Write condensed output file
    taxonomy_coverage_df.write_csv(
        condensed_output_path,
        separator='\t'
    )
    
    logger.info(f"Written genome-wise coverage to {genomewise_output_path}")
    logger.info(f"Written condensed coverage to {condensed_output_path}")


def decompress_genomes_if_needed(genome_metadata_df, temp_dir):
    """
    Decompress genome files if they are gzipped, returning updated metadata.
    
    Args:
        genome_metadata_df: polars DataFrame with genome metadata
        temp_dir: Path object for temporary directory
        
    Returns:
        polars.DataFrame: Updated metadata with decompressed file paths
    """
    logger.info("Checking for compressed genomes...")
    
    decompressed_paths = []
    
    for row in genome_metadata_df.rows(named=True):
        original_path = Path(row['fasta_path'])
        
        if original_path.suffix == '.gz':
            # Create decompressed file in temp directory
            basename = original_path.stem  # Remove .gz
            decompressed_path = temp_dir / basename
            
            logger.debug(f"Decompressing {original_path} to {decompressed_path}")
            
            try:
                with gzip.open(original_path, 'rb') as f_in:
                    with decompressed_path.open('wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                
                decompressed_paths.append(str(decompressed_path))
                
            except Exception as e:
                logger.error(f"Failed to decompress {original_path}: {e}")
                raise
        else:
            # File is already uncompressed, use original path
            decompressed_paths.append(str(original_path))
    
    # Update the dataframe with decompressed paths
    updated_df = genome_metadata_df.with_columns(
        pl.Series('fasta_path', decompressed_paths)
    )
    
    num_decompressed = sum(1 for path in decompressed_paths if temp_dir.as_posix() in path)
    logger.info(f"Decompressed {num_decompressed} genome files")
    
    return updated_df


def execute_simulation_command(command, index):
    """
    Execute a single ART simulation command.
    
    Args:
        command: Shell command to execute
        index: Index for logging purposes
        
    Returns:
        tuple: (index, success, stdout, stderr)
    """
    try:
        output = subprocess.run(
            ["bash", "-o", "pipefail", "-c", command],
            capture_output=True, text=True
        )
        
        if output.returncode != 0:
            return (index, False, output.stdout, output.stderr)
        else:
            return (index, True, output.stdout, output.stderr)
            
    except Exception as e:
        return (index, False, "", str(e))


def concatenate_and_compress_reads(temp_dir, read_pattern, output_path, read_type, threads):
    """
    Concatenate and compress FASTQ files using native tools for speed.
    
    Args:
        temp_dir: Path object for temporary directory
        read_pattern: Glob pattern relative to temp_dir (e.g., 'simulated_reads/*1.fq')
        output_path: Output path for compressed file
        read_type: String identifier ('R1' or 'R2') for logging
        
    Returns:
        tuple: (read_type, success, error_message)
    """
    try:
        # Find all matching read files using pathlib
        read_files = sorted(temp_dir.glob(read_pattern))
        
        if not read_files:
            return (read_type, False, f"No files found matching pattern: {temp_dir / read_pattern}")
        
        logger.debug(f"Found {len(read_files)} files for {read_type} concatenation")

        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Ensure pigz always gets a valid process count.
        pigz_threads = max(1, int(threads) // 2)

        # Use native tools for maximum speed
        files_str = ' '.join(f'"{f}"' for f in read_files)  # Quote filenames for safety
        cmd = f"cat {files_str} | sed 's=/= =' | pigz -p {pigz_threads}  > {output_path}"
        
        logger.debug(f"Concatenation command for {read_type}: {cmd}")
        
        result = subprocess.run(
            ["bash", "-o", "pipefail", "-c", cmd],
            capture_output=True, text=True
        )
        
        if result.returncode != 0:
            return (read_type, False, f"Command failed: {result.stderr}")
        
        return (read_type, True, "")
        
    except Exception as e:
        return (read_type, False, str(e))


def simulate_reads_and_concatenate(simulation_commands, threads, temp_dir, output_read1_path, output_read2_path):
    """
    Execute ART simulations and concatenate results.
    
    Args:
        simulation_commands: List of ART commands to execute
        threads: Number of threads for parallel execution
        temp_dir: Path object for temporary directory
        output_read1_path: Output path for R1 reads
        output_read2_path: Output path for R2 reads
    """
    logger.info(f"Simulating reads for {len(simulation_commands)} genomes...")

    # Snakemake may scale threads down to available cores at runtime.
    threads = max(1, int(threads))
    
    # Execute simulation commands in parallel using ThreadPoolExecutor
    failed_commands = []
    
    with ThreadPoolExecutor(max_workers=threads) as executor:
        # Submit all simulation jobs
        future_to_command = {
            executor.submit(execute_simulation_command, cmd, idx): (idx, cmd)
            for idx, cmd in enumerate(simulation_commands)
        }
        
        # Collect results as they complete
        for future in as_completed(future_to_command):
            idx, cmd = future_to_command[future]
            try:
                index, success, stdout, stderr = future.result()
                
                if not success:
                    logger.error(f"Simulation {index} failed:")
                    logger.error(f"Command: {cmd}")
                    logger.error(f"stdout: {stdout}")
                    logger.error(f"stderr: {stderr}")
                    failed_commands.append((index, cmd, stderr))
                else:
                    logger.debug(f"Simulation {index} completed successfully")
                    
            except Exception as e:
                logger.error(f"Exception in simulation {idx}: {e}")
                failed_commands.append((idx, cmd, str(e)))
    
    # Check if any simulations failed
    if failed_commands:
        logger.error(f"{len(failed_commands)} simulation commands failed:")
        for idx, cmd, error in failed_commands[:3]:  # Show first 3 failures
            logger.error(f"  {idx}: {error}")
        raise Exception(f"Error in {len(failed_commands)} simulation commands.")
        
    logger.info("Read simulation completed successfully")
    
    logger.info("Concatenating and compressing simulated reads...")
    
    with ThreadPoolExecutor(max_workers=2) as executor:
        # Submit both concatenation jobs
        future_r1 = executor.submit(
            concatenate_and_compress_reads, 
            temp_dir,
            'simulated_reads/*1.fq', 
            output_read1_path, 
            'R1',
            threads
        )
        future_r2 = executor.submit(
            concatenate_and_compress_reads, 
            temp_dir,
            'simulated_reads/*2.fq', 
            output_read2_path, 
            'R2',
            threads
        )
        
        # Collect results
        failed_concatenations = []
        
        for future, read_type in [(future_r1, "R1"), (future_r2, "R2")]:
            try:
                read_type_result, success, error_message = future.result()
                
                if not success:
                    logger.error(f"{read_type} concatenation failed: {error_message}")
                    failed_concatenations.append((read_type, error_message))
                else:
                    logger.info(f"{read_type} concatenation and compression completed successfully")
                    
            except Exception as e:
                logger.error(f"Exception in {read_type} concatenation: {e}")
                failed_concatenations.append((read_type, str(e)))
    
    # Check if any concatenations failed
    if failed_concatenations:
        error_msg = "; ".join([f"{read_type}: {error}" for read_type, error in failed_concatenations])
        raise Exception(f"Error in concatenation: {error_msg}")
        
    logger.info(f"Concatenated reads written to {output_read1_path} and {output_read2_path}")



def generate_mock_community(coverage_file_path, genome_list_path, metadata_file_path, 
                          art_binary_path, threads, output_read1_path, output_read2_path,
                          condensed_output_path, genomewise_output_path):
    """
    Generate a mock microbial community with simulated reads.
    
    Args:
        coverage_file_path: Path to coverage values file
        genome_list_path: Path to genome list file
        metadata_file_path: Path to metadata file
        art_binary_path: Path to ART binary
        threads: Number of threads to use
        output_read1_path: Output path for R1 reads
        output_read2_path: Output path for R2 reads
        condensed_output_path: Output path for condensed coverage
        genomewise_output_path: Output path for genome-wise coverage
    """
    logger.info("Starting mock community generation...")
    
    # Load input data
    coverage_values = load_coverage_data(coverage_file_path)
    genome_metadata_df = load_genome_metadata(genome_list_path, metadata_file_path)
    
    # Convert output paths to absolute paths before entering temp directory context
    output_read1_path = Path(output_read1_path).resolve()
    output_read2_path = Path(output_read2_path).resolve()
    condensed_output_path = Path(condensed_output_path).resolve()
    genomewise_output_path = Path(genomewise_output_path).resolve()
    
    # Use temporary directory for simulation
    with tempfile.TemporaryDirectory() as temp_dir_str:
        temp_dir = Path(temp_dir_str)
        
        # Create simulated_reads subdirectory
        simulated_reads_dir = temp_dir / 'simulated_reads'
        simulated_reads_dir.mkdir(exist_ok=True)
        logger.debug(f"Working with temporary directory: {temp_dir}")
        
        # Decompress genomes if necessary
        decompressed_genome_metadata = decompress_genomes_if_needed(genome_metadata_df, temp_dir)
        
        # Generate simulation commands with decompressed paths
        simulation_commands, genome_metadata_with_coverage = generate_simulation_commands(
            decompressed_genome_metadata, coverage_values, art_binary_path, 
            READ_LENGTH, FRAGMENT_LENGTH, FRAGMENT_STDEV, temp_dir
        )
        
        # Write output files (use original metadata for paths)
        coverage_file_basename = Path(coverage_file_path).name
        write_output_files(
            genome_metadata_with_coverage, coverage_file_basename,
            condensed_output_path, genomewise_output_path
        )
        
        # Simulate reads and concatenate
        simulate_reads_and_concatenate(
            simulation_commands, threads, temp_dir, output_read1_path, output_read2_path
        )
    
    logger.info("Mock community generation completed successfully")

# main
# ----
def main():
    """Main function to orchestrate mock community generation."""
    generate_mock_community(
        COVERAGE_FILE, GENOME_LIST, METADATA_FILE, ART_BINARY, THREADS,
        OUTPUT_READ1, OUTPUT_READ2, OUTPUT_CONDENSED, OUTPUT_GENOMEWISE_COVERAGE
    )

# run
# ---
if __name__ == "__main__":
    main()
    