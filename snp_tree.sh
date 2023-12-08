## SNP tree building
## Author: Jie.Lu@dshs.texas.gov
## Version: 1.0
## Create a directory naming after the request date and copy the corresponding assembled contigs to the dir
## Usage: snp_tree.sh <speciesDir/DateOfRequest>


baseDir="/home/dnalab/cluster_analysis"
fastaDir=$baseDir/$1
outDir=$baseDir/$1/result
cd $outDir
nextflow run /pipeline/ARLN/Grandeur/grandeur.nf -c /pipeline/ARLN/Grandeur/configs/tree.config -profile singularity --fastas $fastaDir --outdir $outDir
rm -r $ourDir/work
