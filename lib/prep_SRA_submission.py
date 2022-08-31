#!/usr/bin/env python3

import pandas as pd
from glob import glob
from os import path
from datetime import date

def prep_SRA_submission(qc_results, run_name):
    reads_dir = "/home/dnalab/reads/{}/".format(run_name)
    metadata = pd.read_csv("/home/jiel/bin/armadillo/template/SRA_metadata_template.txt", sep="\t", header=0, index_col=None)
    attribute = pd.read_csv("/home/jiel/bin/armadillo/template/attribute_template.txt", sep="\t", header=0, index_col=None)
    results = pd.read_csv(qc_results, sep="\t", header=0, index_col=None)

    for i, row in results.iterrows():
        if row["Status"] == "Complete":
            sample_id = row["sample_id"]
            fastq_files = []
            fastqs = reads_dir + sample_id + "*"
            for fastq in glob(fastqs):
                fastq_files.append(path.basename(fastq))
    
            new_row_metadata = {"sample_name": sample_id, "library_ID": sample_id, "title": "Illumina sequencing of {}".format(sample_id), 
               "library_strategy": "WGS", "library_source": "GENOMIC",	"library_selection": "RANDOM", "library_layout": "PAIRED",	
               "platform": "ILLUMINA",	"instrument_model": "Illumina MiSeq", "design_description": "Illumina DNA Prep", "filetype": "FASTQ",
               "filename": fastq_files[0],	"filename2": fastq_files[1], "filename3": "", "filename4": "",	"assembly": "",	"fasta_file": ""}
            metadata = metadata.append(new_row_metadata, ignore_index = True)

            new_row_attr = {"*sample_name": sample_id, "bioproject_accession": "PRJNA288601", "*organism": row["kraken2_top_species"], 
                            "*collected_by": "Missing", "*collection_date": date.today().year, "*geo_loc_name": "USA", "*host":"Homo sapiens",	
                            "*host_disease":"missing",	"*isolate": "whole organism", "*isolation_source":"missing", "*lat_lon": "missing"}
            attribute = attribute.append(new_row_attr, ignore_index = True)
    metadata.to_csv(run_name + "_SRA_metadata.tsv", sep = "\t", index = False)
    attribute.to_csv(run_name + "_SRA_attribute.tsv", sep = "\t", index = False)

if __name__ == "__main__":
    prep_SRA_submission("qc_results.tsv", "AR_220819_M03431")