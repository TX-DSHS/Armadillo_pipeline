#!/usr/bin/env python3
# Armadillo pipeline 
# Author: Jie.Lu@dshs.texas.gov
version = 2.0-10/3/2023

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
from lib.prep_SRA_submission_phoenix import prep_SRA_submission
from lib.presence_gene_family import check_gene_family

# def updateSQLiteTable(results):
#     conn = sqlite3.connect("/home/dnalab/database/arln.sqlite")
#     results.to_sql("results", conn, if_exists="append")
#     # cur = conn.cursor()
#     # cur.execute("SELECT * FROM results")
#     # rows = cur.fetchall()
#     # for row in rows:
#     #     print(row)
#     conn.close()
#     return

# def readMetadata():
#     # Read sqlite query results into a pandas DataFrame
#     conn = sqlite3.connect("/home/dnalab/database/arln.sqlite")
#     df = pd.read_sql_query("SELECT * from metadata", conn)

#     # Verify that result of SQL query is stored in the dataframe
#     #print(df.head())
#     conn.close()
#     return df

my_parser = argparse.ArgumentParser()
my_parser.add_argument("-i", help = 'Grandeur output: grandeur_results.tsv', default = "results/Phoenix_Summary.tsv")
my_parser.add_argument("-r", help = 'Run name')

args = my_parser.parse_args()
phoenixSummaryFile = args.i
run_name = args.r

logging.basicConfig(filename = 'armadillo.log', filemode = 'a', level = logging.DEBUG)
logging.info('Armadillo {} starting on run {}'.format(version, run_name))
logging.info(str(date.today()))

#####################################################################
# Read Phoenix_summary.tsv and extract Hai-seq ID from sample name
#######################################################################
results = pd.read_csv(phoenixSummaryFile, sep="\t", header=0, index_col=None)
results["HAIseq_ID"] = results["ID"].str[:12]
results["run_name"] = run_name
results['GAMMA_Beta_Lactam_Resistance_Genes'] = results['GAMMA_Beta_Lactam_Resistance_Genes'].fillna("Not detected")
results['GAMMA_Other_AR_Genes'] = results['GAMMA_Other_AR_Genes'].fillna("Not detected")


# # update the SQLite database
# try:
#     updateSQLiteTable(results)
# except:
#     pass

########################################################
# Check the presence of non-OXA family genes
#########################################################
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

for gene_list in results["GAMMA_Beta_Lactam_Resistance_Genes"]:
    for carb_gene in carb_genes:
        checklist[carb_gene] = "NOT_DETECTED"

    gene_list = gene_list.replace('"', '').split(",")
    gene_list = [ i.split("_")[0] for i in gene_list]

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

results.sort_values(by="HAIseq_ID", ascending=True, inplace=True)

#####################################################################
# Generate SRA submission files
#####################################################################

results_to_sra = results[results["Auto_QC_Outcome"] == "PASS"]
results_to_sra = results_to_sra[results_to_sra["ID"].str.contains('CON') == False] 
results_metadata = prep_SRA_submission(results_to_sra, run_name)
print(results_metadata)

passSamples = list(results_metadata[results_metadata["Auto_QC_Outcome"] == "PASS"]["ID"])
passSample_Ids = list(results_metadata[results_metadata["Auto_QC_Outcome"] == "PASS"]["HAIseq_ID"])
passSample_name = list(results_metadata[results_metadata["Auto_QC_Outcome"] == "PASS"]["Species"])
passSample_mlst = list(results_metadata[results_metadata["Auto_QC_Outcome"] == "PASS"]["MLST_1"])
passSample_specimen_id = list(results_metadata[results_metadata["Auto_QC_Outcome"] == "PASS"]["KEY"])
passSample_beta_lactam = list(results_metadata[results_metadata["Auto_QC_Outcome"] == "PASS"]["GAMMA_Beta_Lactam_Resistance_Genes"])
passSample_other_AR = list(results_metadata[results_metadata["Auto_QC_Outcome"] == "PASS"]["GAMMA_Other_AR_Genes"])

#amrheader = ["Gene symbol", "Sequence name"]

data_uri = base64.b64encode(open('/home/jiel/bin/armadillo/DSHS_Banner.png', 'rb').read()).decode('utf-8')
img_tag = '<img src="data:image/png;base64,{0}">'.format(data_uri)
footnote = open("/home/jiel/bin/armadillo/footnote.txt", 'r')
footnote = footnote.readlines()
footnote = ("\n").join(footnote)

reads_dir = "/home/dnalab/reads/{}/".format(run_name)
subprocess.run(["mkdir", "-p", "SRA_seq"])

refseq = pd.read_csv("/home/dnalab/ReferenceGeneCatalog_3.11_20230417.txt", sep="\t", header=0, index_col=None)
refseq['allele'] = refseq['allele'].fillna(refseq['gene_family'])



for s, id, name, mlst, specimen_id, beta_lactam, other_AR in zip(passSamples, passSample_Ids, passSample_name, passSample_mlst, passSample_specimen_id, passSample_beta_lactam, passSample_other_AR):

    # create pdf reports for each passed sample    

    # check the presence of OXA family genes

    outname = id +"_amrfinder_plus_report.txt"
    if beta_lactam != "Not detected" and other_AR != "Not detected":
        all_amr_genes = beta_lactam.split(",") + other_AR.split(",")
        

    elif beta_lactam != "Not detected" and other_AR == "Not detected":
        all_amr_genes = beta_lactam.split(",")
        
    
    elif beta_lactam == "Not detected" and other_AR != "Not detected":
        all_amr_genes = other_AR.split(",")
    
    else:
        all_amr_genes = ["Not detected"]
        
    amrheader = ["Gene_symbol", "product_name"]
    all_amr_genes = [ i.split("_")[0] for i in all_amr_genes]
    all_amr_genes = pd.DataFrame({"Gene_symbol": all_amr_genes})
    all_amr_genes = pd.merge(all_amr_genes, refseq, left_on = "Gene_symbol", right_on = "allele", how = "left")
    all_amr_genes["product_name"] = all_amr_genes["product_name"].fillna("Description not found")
    all_amr_genes = all_amr_genes[["Gene_symbol", "product_name"]]
    all_amr_genes = all_amr_genes.drop_duplicates()
    
    all_amr_genes.sort_values(by="Gene_symbol", ascending=True, inplace=True)
    all_amr_genes.to_csv(outname, sep = "\t", columns = amrheader, index = False)

    try:
        oxa_families = check_gene_family(all_amr_genes["product_name"])
        for family in oxa_families:
            results.loc[results["HAIseq_ID"] == id, family] = oxa_families[family]
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
    
    # link fastq files to SRA_seq folder
    fastqs = reads_dir + id + "*"
    for fastq in glob(fastqs):
        #print(fastq)
        filelink = "SRA_seq/"+ path.basename(fastq)
        subprocess.run(["ln", "-s", fastq, filelink])


    #copy contigs fastas to cluster folder on S3
    print(specimen_id)
    genus = name.split(" ")[0]
    species = name.split(" ")[1]
    contig_fasta = "results/" + s + "/assembly"+ s + ".contigs.fa.gz"
    fasta_path = "/home/dnalab/cluster/" + genus + "_" + species
    fasta_name = specimen_id + '_' + s + "_contigs.fa.gz"
    print(contig_fasta, fasta_name)
    # if not path.exists(fasta_path):
    #    system("mkdir {}".format(fasta_path))
    # else:
    #    system("cp {} {}/{}".format(contig_fasta, fasta_path, fasta_name))
    #    system("aws s3 cp {} s3://804609861260-bioinformatics-infectious-disease/cluster/{}_{}/{} --region us-gov-west-1".format(contig_fasta, genus, species, fasta_name)) 

#####################################
# write results to qc_results.xlsx
#####################################
# Header of Phoenxi_summary.tsv
# HAIseq_ID	Auto_QC_Outcome	Warning_Count	Estimated_Coverage	Genome_Length	Assembly_Ratio_(STDev)	#_of_Scaffolds_>500bp	GC_%	
# Species	Taxa_Confidence	Taxa_Coverage	Taxa_Source	Kraken2_Trimd	Kraken2_Weighted	MLST_Scheme_1	MLST_1	MLST_Scheme_2	
# MLST_2	GAMMA_Beta_Lactam_Resistance_Genes	GAMMA_Other_AR_Genes	AMRFinder_Point_Mutations	Hypervirulence_Genes	
# Plasmid_Incompatibility_Replicons	Auto_QC_Failure_Reason

column = ["HAIseq_ID", "run_name", "Species", "Auto_QC_Outcome", "Auto_QC_Failure_Reason", 
          "Estimated_Coverage", "Genome_Length", "Assembly_Ratio_(STDev)", "#_of_Scaffolds_>500bp", "GC_%", 
          "blaKPC", "blaNDM", "blaOXA-48", "blaVIM", "blaIMP", "blaOXA-23", "blaOXA-24/40", "blaOXA-58", 
          "Hypervirulence_Genes", "MLST_1"
         ]
results.to_csv("qc_results.tsv", sep = "\t", columns = column, index = False)
results.to_excel("qc_results.xlsx", header = True, columns = column, index = False)
logging.info('Armadillo finished')