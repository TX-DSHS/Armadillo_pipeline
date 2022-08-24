#!/usr/bin/env python3

import sys
import re
import pandas as pd
import numpy as np
import pdfkit
from prettytable import PrettyTable
print("running armadillo...")
grandeurSummaryFile = sys.argv[1]
results = pd.read_csv(grandeurSummaryFile, sep="\t", header=0, index_col=None)
passQC = np.where((results["cg_coverage"] >= 40) & (results["quast_contigs"] < 200), "Complete", "Failed_QC")
QCtag1 = np.where(results["cg_coverage"] < 40, "Low_Coverage(<40)", "")
QCtag2 = np.where(results["quast_contigs"] >=200, "Too_many_contigs(>=200)", "")
QCtag = []
for a,b in zip(QCtag1, QCtag2):
    if a == "" and b == "":
        QCtag.append("")
    else:
        QCtag.append(a+";"+b)
results["Status"] = passQC
results["QCtag"] = QCtag
#print(results)
#print(QCtag)
#sample_id       sample  seqyclean_pairs_kept    seqyclean_percent_kept  fastqc_1_reads  fastqc_2_reads  mash_genome_size        mash_coverage   mash_genus      mash_species    mash_full       mash_pvalue
#     mash_distance   fastani_ref_top_hit     fastani_ani_score       fastani_per_aligned_seq_matches quast_gc_%      quast_contigs   quast_N50       quast_length    cg_average_read_length  cg_average_quality      cg_coverage     ref_genome_length       amr_genes       virulence_genes
checklist = {}
carb_genes = ["blaKPC", "blaNDM", "blaOXA-48", "blaVIM", "blaIMP", "blaOXA-23", "blaOXA-24", "blaOXA-40", "blaOXA-24/40", "blaOXA-58"]
carb_gene_labels = {}
for carb_gene in carb_genes:
    carb_gene_labels[carb_gene] = []

for gene_list in results["amr_genes"]:
    for carb_gene in carb_genes:
        checklist[carb_gene] = "No"
    print(gene_list)
    gene_list = gene_list.replace('"', '').split(",")
    for gene in gene_list:
        print(gene)
        if gene in checklist:
            checklist[gene] = "Yes"
        elif gene.startswith("blaKPC"):
            checklist["blaKPC"] = "Yes"
        elif gene.startswith("blaNDM"):
            checklist["blaNDM"] = "Yes"
        elif gene.startswith("blaVIM"):
            checklist["blaVIM"] = "Yes"
        elif gene.startswith("blaIMP"):
            checklist["blaIMP"] = "Yes"
        
    print(checklist)
    for gene in checklist:
        carb_gene_labels[gene].append(checklist[gene])

print(carb_gene_labels)
for gene in carb_gene_labels:
    results[gene] = carb_gene_labels[gene]

header = ["sample_id", "sample", "mash_genus", "mash_species", "Status", "QCtag", "cg_coverage", "quast_contigs", "blaKPC", "blaNDM", "blaOXA-48", "blaVIM", "blaIMP", "blaOXA-23", "blaOXA-24", "blaOXA-40", "blaOXA-24/40", "blaOXA-58"]
results.to_csv("qc_results.tsv", sep = "\t", columns = header, index = False)

passSamples = list(results[results["Status"] == "Complete"]["sample"])
amrheader = ["Gene symbol", "Sequence name"]
print(passSamples)
for s in passSamples:
    amrFile = "ncbi-AMRFinderplus/"+ s + "_amrfinder_plus.txt"
    df = pd.read_csv(amrFile, sep="\t", header=0, index_col=None)
    outname = s+"_amrfinder_plus_report.txt"
    #Name	Protein identifier	Contig id	Start	Stop	Strand	Gene symbol	Sequence name	Scope	Element type	Element subtype	Class	Subclass	Method	Target length	Reference sequence length	% Coverage of reference sequence	% Identity to reference sequence	Alignment length	Accession of closest sequence	Name of closest sequence	HMM id	HMM description
    df.to_csv(outname, sep = "\t", columns = amrheader, index = False)
    # open csv file
    a = open(outname, 'r')
    # read the csv file
    a = a.readlines()
    # Separating the Headers
    l1 = a[0]
    l1 = l1.split('\t')

    # headers for table
    t = PrettyTable(l1)

    # Adding the data
    for i in range(1, len(a)) :
        t.add_row(a[i].split('\t'))

    code = t.get_html_string()
    html_file = open(s+"_amrfinder_plus_report.html", 'w')
    html_file = html_file.write(code)
    pdfkit.from_file(s+"_amrfinder_plus_report.html", s+"_amrfinder_plus_report.pdf")

print("armadillo completed")