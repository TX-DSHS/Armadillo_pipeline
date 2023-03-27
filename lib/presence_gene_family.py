#!/usr/bin/env python3
# Armadillo pipeline 
# Author: Jie.Lu@dshs.texas.gov

import re

def check_gene_family(gene_list):
    oxa_families = {"blaOXA-23":"NOT_DETECTED",
                    "blaOXA-24/40": "NOT_DETECTED", 
                    "blaOXA-58": "NOT_DETECTED",
                    "blaOXA-48": "NOT_DETECTED"}
    for amrgene in gene_list:
        if re.search(r"\bOXA-24\b", amrgene) or re.search(r"\bOXA-40\b", amrgene):
            oxa_families["blaOXA-24/40"] = "DETECTED"
        elif re.search(r"\bOXA-58\b", amrgene):
            oxa_families["blaOXA-58"] = "DETECTED"
        elif re.search(r"\bOXA-23\b", amrgene):
            oxa_families["blaOXA-23"] = "DETECTED"
        elif re.search(r"\bOXA-48\b", amrgene):
            oxa_families["blaOXA-48"] = "DETECTED"
    return oxa_families

if __name__ == "__main__":
    print(check_gene_family(["OXA-241 gene"]))
    print(check_gene_family(["OXA-24 gene", "OXA-23 gene", "family OXA-58 gene"]))
    print(check_gene_family(["OXA-24", "OXA-40"]))
    print(check_gene_family(["family OXA-4011", "family OXA-58"]))
    print(check_gene_family(["OXA-400", "OXA-581"]))
        