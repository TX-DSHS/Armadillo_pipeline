#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np
from os import path, system

def id_match(qc_result, metadata):
    # check if kraken2 ID matches workbook ID. Assuming metadata file has been validated.
    results = pd.read_csv(qc_result, sep="\t", header=0, index_col=None)
    metadata = pd.read_csv(metadata, sep="\t", header=0, index_col=None)
    results = pd.merge(results, metadata, left_on = "sample_id", right_on = "Sample_ID", how = "left")
    results["ID_Matched"] = np.where(results["kraken2_top_species"] == results["Species"], "Matched", "Warning: IDs do not match")
    
    return results

if __name__ == "__main__":
    print(id_match("qc_results.tsv", "dummy_metadata.tsv"))
