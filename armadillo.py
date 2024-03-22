#!/usr/bin/env python3
# Armadillo pipeline 
# Author: Jie.Lu@dshs.texas.gov
version = 0.1-2/21/2023

import sys
import argparse
import re
import pandas as pd
import numpy as np
import sqlite3
from datetime import date
import pdfkit
from prettytable import PrettyTable
import base64
from os import path, system
import logging
import subprocess
from glob import glob
from lib.prep_SRA_submission import prep_SRA_submission
from lib.presence_gene_family import check_gene_family

def updateSQLiteTable(results):
    conn = sqlite3.connect("/home/dnalab/database/arln.sqlite")
    results.to_sql("results", conn, if_exists="append")
    # cur = conn.cursor()
    # cur.execute("SELECT * FROM results")
    # rows = cur.fetchall()
    # for row in rows:
    #     print(row)
    conn.close()
    return

def readMetadata():
    # Read sqlite query results into a pandas DataFrame
    conn = sqlite3.connect("/home/dnalab/database/arln.sqlite")
    df = pd.read_sql_query("SELECT * from metadata", conn)

    # Verify that result of SQL query is stored in the dataframe
    #print(df.head())
    conn.close()
    return df

my_parser = argparse.ArgumentParser()
my_parser.add_argument("-i", help = 'Grandeur output: grandeur_results.tsv', default = "grandeur_results.tsv")
my_parser.add_argument("-r", help = 'Run name')
my_parser.add_argument("-c", help = 'Minimal coverage (defaul = 30)', default = 30)
args = my_parser.parse_args()
grandeurSummaryFile = args.i
run_name = args.r
min_coverage = args.c

logging.basicConfig(filename = 'armadillo.log', filemode = 'a', level = logging.DEBUG)
logging.info('Armadillo v0.1 starting on run {}'.format(run_name))
logging.info(str(date.today()))

results = pd.read_csv(grandeurSummaryFile, sep="\t", header=0, index_col=None)

# Check coverage and number of contigs
passQC = np.where((results["cg_coverage"] >= min_coverage), "Complete", "Failed_QC")
QCtag1 = np.where(results["cg_coverage"] < min_coverage, "Low_Coverage(<{})".format(min_coverage), "")
QCtag2 = np.where(results["quast_contigs"] >=200, "Warining:large_num_contigs(>=200)", "")
QCtag = []
for a,b in zip(QCtag1, QCtag2):
    if a == "" and b == "":
        QCtag.append("")
    elif a == "" and b != "":
        QCtag.append(b)
    elif a != "" and b == "":
        QCtag.append(a)
    else:
        QCtag.append(a+";"+b)
results["Status"] = passQC
results["QCtag"] = QCtag
results["analysis_date"] = [str(date.today())] * len(results)
results["run_name"] = [run_name] * len(results)

# update the SQLite database
try:
    updateSQLiteTable(results)
except:
    pass

# Check genome size and GC content
ncbi_genome_size = pd.read_csv("/home/jiel/bin/armadillo/genome_size.txt", sep = "\t", header = 0)
results["MASH_ID"] = results["mash_genus"] + " " + results["mash_species"]
results = pd.merge(results, ncbi_genome_size, left_on = "kraken2_top_species", right_on = "Organism_ID", how = "left")

#sample_id       sample  seqyclean_pairs_kept    seqyclean_percent_kept  fastqc_1_reads  fastqc_2_reads  mash_genome_size        mash_coverage   mash_genus      mash_species    mash_full       mash_pvalue
#     mash_distance   fastani_ref_top_hit     fastani_ani_score       fastani_per_aligned_seq_matches quast_gc_%      quast_contigs   quast_N50       quast_length    cg_average_read_length  cg_average_quality      cg_coverage     ref_genome_length       amr_genes       virulence_genes   mlst
checklist = {}
carb_genes = ["blaKPC", "blaNDM", "blaOXA-48", "blaVIM", "blaIMP", "blaOXA-23", "blaOXA-24/40", "blaOXA-58"]
oxa_family = {}
with open("/home/jiel/bin/armadillo/oxa_family.tsv", "r") as f:
    for line in f:
        line = line.rstrip()
        family, oxa_gene = line.split("\t")
        oxa_family[oxa_gene] = family

carb_gene_labels = {}
for carb_gene in carb_genes:
    carb_gene_labels[carb_gene] = []

for gene_list in results["amr_genes"]:
    for carb_gene in carb_genes:
        checklist[carb_gene] = "NOT_DETECTED"

    gene_list = gene_list.replace('"', '').split(",")
    for gene in gene_list:
        if gene in oxa_family:
            checklist[oxa_family[gene]] = "DETECTED"
        elif gene.startswith("blaKPC"):
            checklist["blaKPC"] = "DETECTED"
        elif gene.startswith("blaNDM"):
            checklist["blaNDM"] = "DETECTED"
        elif gene.startswith("blaVIM"):
            checklist["blaVIM"] = "DETECTED"
        elif gene.startswith("blaIMP"):
            checklist["blaIMP"] = "DETECTED"


    #print(checklist)
    for gene in checklist:
        carb_gene_labels[gene].append(checklist[gene])

#print(carb_gene_labels)
for gene in carb_gene_labels:
    results[gene] = carb_gene_labels[gene]

column = ["sample_id", "run_name", "kraken2_top_species", "MASH_ID", "Status", "QCtag", "cg_coverage", 
          "quast_contigs", "quast_length", "quast_gc_%", "blaKPC", "blaNDM", "blaOXA-48", "blaVIM", 
          "blaIMP", "blaOXA-23", "blaOXA-24/40", "blaOXA-58", "virulence_genes", "mlst"]

results.sort_values(by="sample_id", ascending=True, inplace=True)

# Generate SRA submission files
results_to_sra = results[results["Status"] == "Complete"]
results_to_sra = results_to_sra[results_to_sra["sample_id"].str.contains('CON') == False] 
results_metadata = prep_SRA_submission(results_to_sra, run_name)

passSamples = list(results_metadata[results_metadata["Status"] == "Complete"]["sample"])
passSample_Ids = list(results_metadata[results_metadata["Status"] == "Complete"]["sample_id"])
passSample_name = list(results_metadata[results_metadata["Status"] == "Complete"]["kraken2_top_species"])
passSample_mlst = list(results_metadata[results_metadata["Status"] == "Complete"]["mlst"])
passSample_specimen_id = list(results_metadata[results_metadata["Status"] == "Complete"]["KEY"])
amrheader = ["Gene symbol", "Sequence name"]

data_uri = base64.b64encode(open('/home/jiel/bin/armadillo/DSHS_Banner.png', 'rb').read()).decode('utf-8')
img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
footnote = open("/home/jiel/bin/armadillo/footnote.txt", 'r')
footnote = footnote.readlines()
footnote = ("\n").join(footnote)

reads_dir = "/home/dnalab/reads/{}/".format(run_name)
subprocess.run(["mkdir", "-p", "SRA_seq"])

for s, id, name, mlst, specimen_id in zip(passSamples, passSample_Ids, passSample_name, passSample_mlst, passSample_specimen_id):
    # link fastq files to SRA_seq folder
    fastqs = reads_dir + s + "*"
    for fastq in glob(fastqs):
        #print(fastq)
        filelink = "SRA_seq/"+ path.basename(fastq)
        subprocess.run(["ln", "-s", fastq, filelink])

    # create pdf reports for each passed sample    
    genus = name.split(" ")[0]
    species = name.split(" ")[1]
    amrFile = "ncbi-AMRFinderplus/"+ s + "_amrfinder_plus.txt"
    df = pd.read_csv(amrFile, sep="\t", header=0, index_col=None)
    outname = id +"_amrfinder_plus_report.txt"
    #Name	Protein identifier	Contig id	Start	Stop	Strand	Gene symbol	Sequence name	Scope	Element type	Element subtype	Class	Subclass	Method	Target length	Reference sequence length	% Coverage of reference sequence	% Identity to reference sequence	Alignment length	Accession of closest sequence	Name of closest sequence	HMM id	HMM description
    df.to_csv(outname, sep = "\t", columns = amrheader, index = False)
    try:
        oxa_families = check_gene_family(df["Name of closest sequence"])
        for family in oxa_families:
            results.loc[results["sample_id"] == id, family] = oxa_families[family]
    except:
        pass
    # open csv file
    a = open(outname, 'r')
    # read the csv file
    a = a.readlines()
    # Separating the Headers
    l1 = a[0]
    l1 = l1.split('\t')

    # headers for table
    t = PrettyTable(border=True, header=True, padding_width=5)
    t.field_names=l1
    # Adding the data
    for i in range(1, len(a)) :
        t.add_row(a[i].split('\t'))
    
    code = t.get_html_string(attributes={'border': 1, 'style': 'border-width: 1px; border-collapse: collapse; font-size:25px'})
    html_file = open(id + "_amrfinder_plus_report.html", 'w')
    timestamp = str(date.today())
    html_file.write("<h1 style=\"text-align:center\">Next-Generation Sequencing (NGS) Analysis Report</h1>\n")
    html_file.write(img_tag)
    html_file.write("<h2>Report date: " + timestamp + "</h2>\n")
    html_file.write("<h2>Sample ID: " + id + "</h2>\n")
    #html_file.write("<h2>TX ID: " + specimen_id + "</h2>\n")
    html_file.write("<h2>MLST: " + mlst + "</h2>\n")
    html_file.write("<h3>Table: Antimicrobial Resistance genes identified</h3>\n")
    html_file.write(code)
    html_file.write("<br><footer>"+footnote+"</footer>")
    
    options = {
    'page-size': 'Letter',
    'margin-top': '0.75in',
    'margin-right': '0.75in',
    'margin-bottom': '0.75in',
    'margin-left': '0.75in',
    'encoding': "UTF-8",
    'footer-right': '[page] of [topage]',
    'no-outline': None
}
    html_file.close()
    pdfkit.from_file(id +"_amrfinder_plus_report.html", id + "_amrfinder_plus_report.pdf", options = options)

    #copy contigs fastas to cluster
    print(specimen_id)
    contig_fasta = "contigs/"+ s + "_contigs.fa"
    fasta_path = "/home/dnalab/cluster/" + genus + "_" + species
    fasta_name = specimen_id + '_' + s + "_contigs.fa"
    if not path.exists(fasta_path):
       system("mkdir {}".format(fasta_path))
    else:
       system("cp {} {}/{}".format(contig_fasta, fasta_path, fasta_name))
       system("aws s3 cp {} s3://804609861260-bioinformatics-infectious-disease/cluster/{}_{}/{} --region us-gov-west-1".format(contig_fasta, genus, species, fasta_name)) 

results.to_csv("qc_results.tsv", sep = "\t", columns = column, index = False)
results.to_excel("qc_results.xlsx", header = True, columns = column, index = False)
logging.info('Armadillo finished')