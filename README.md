To run these benchmarks, first setup the master environment:
```
mamba env create -f env.yml -p ./env
conda activate ./env
```

Update this section of the config with your own paths.
`taxonomy_map` should be a .tsv file with the columns 'acc' (the genome accession or identifier string) and 'parsed_taxon_string' (the taxonomy string used to generate the metapackage)

```
# input paths
genome_fasta_glob_path: '/path/to/reference/database/*.fasta.gz'
singlem_metapackage: "/path/to/metapackage.smpkg/"
taxonomy_map: '/path/to/taxonomy_mapping.tsv'
```

Then run the benchmarks with different numbers of threads:
```
snakemake -c10 --use-conda
```
