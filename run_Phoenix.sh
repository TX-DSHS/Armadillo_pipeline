# This script is used to run the PhoeNix pipeline for the ARLN project.
# The script will download the fastq files from s3, run the PhoeNix pipeline, and then run the Armadillo pipeline.
# The results will be zipped and uploaded to s3.
# Usage: bash run_Phoenix.sh <run_id>
# Example: bash run_Phoenix.sh AR_240321_M04922
# Note: The run_id should be the same as the name of the zipped fastq files in s3.
# The fastq files should be zipped and uploaded to s3 in the following format: s3://804609861260-bioinformatics-infectious-disease/ARLN/RAW_RUNS/<run_id>.zip
# The results will be uploaded to s3 in the following format: s3://804609861260-bioinformatics-infectious-disease/ARLN/ANALYSIS_RESULTS/<run_id>_result.zip
# The logs and reports will be uploaded to s3 in the following format: s3://804609861260-bioinformatics-infectious-disease/ARLN/REPORT/<run_id>_report.zip
# The script will also check the size of the fastq files, if the size is < 1Mb, the files will be moved to a folder named "small_size_fastq".
# Author: Jie.Lu@dshs.texas.gov
# All rights reserved.

# If no argument is provided, print the usage
if [ $# -eq 0 ]; then
  echo "Usage: bash run_Phoenix.sh <run_id>"
  exit 1
fi

version="2.1-04/02/2024"
# Read the aws bucket name from file aws_bucket.txt
aws_bucket=$(cat aws_bucket.txt)
#aws_bucket="s3://804609861260-bioinformatics-infectious-disease"
basedir=$PWD
refGenCatlog="/ReferenceGeneCatalog_3.12_20240205.txt"
phoenix_version="v2.1.1"

mkdir -p $basedir/results
mkdir -p $basedir/reads

rm -rf $basedir/results/$1
mkdir -p $basedir/results/$1
echo "Running Armadillo version: "$version > $basedir/results/$1/armadillo.log

# # Copy and unzip the fastq files from s3
aws s3 cp $aws_bucket/ARLN/RAW_RUNS/$1.zip $basedir/reads/zip --region us-gov-west-1
# if upload failed, upload the log file to aws, exit

if [ $? -ne 0 ]; then
  echo "Failed to download the fastq files from s3." >> $basedir/results/$1/armadillo.log
  exit 1
fi

rm -rf $basedir/reads/$1
mkdir $basedir/reads/$1
unzip -j $basedir/reads/zip/$1.zip -d $basedir/reads/$1

# Check if the file size is < 1Mb, if yes then move to a folder
mkdir $basedir/reads/$1/small_size_fastq
echo "Checking the size of the file..." >> $basedir/results/$1/armadillo.log
touch $basedir/results/$1/$1_failed_file_size.log
for fastq in $basedir/reads/$1/*.gz; do
  myfilesize=$(stat --format=%s $fastq)
  if [ $myfilesize -lt 1000000 ]; then
    mv $fastq $basedir/reads/$1/small_size_fastq
    echo $fastq  >> $basedir/results/$1/failed_file_size.log
  fi
done

# # If all the files are smaller than 1Mb, then exit
# if [ ! -f $basedir/reads/$1/*.gz ]; then
#   echo "All the files are smaller than 1Mb, please check the fastq files." >> $basedir/results/$1/armadillo.log
#   exit 1
# fi

# Run PhoeNix pipeline
reads_path=$basedir/reads/$1
out_path=$basedir/results/$1
$basedir/phoenix/bin/create_samplesheet.sh $reads_path > $reads_path/samplesheet.csv
source $basedir/miniconda3/etc/profile.d/conda.sh
conda activate nextflow
cd $out_path
mkdir -p $basedir/singularity/Phoenix
export NXF_SINGULARITY_CACHEDIR=$basedir/singularity/Phoenix
nextflow run cdcgov/phoenix -r $phoenix_version -profile singularity -entry PHOENIX --input $reads_path/samplesheet.csv --kraken2db $basedir/k2_standard_08gb_20230605 --output $out_path
conda deactivate
# if nextflow failed, exit
if [ $? -ne 0 ]; then
  echo "Failed to run the PhoeNix pipeline." >> $basedir/results/$1/armadillo.log
  exit 1
fi

rm -r $out_path/work
rm $basedir/results/zip/$1.zip

# Run Armadillo pipeline
python3 $basedir/armadillo_phoenix.py -r $1 -d $basedir -a $aws_bucket -c $refGenCatlog >> $basedir/results/$1/armadillo.log

# if armadillo failed, exit
if [ $? -ne 0 ]; then
  echo "Failed to run the Armadillo pipeline." >> $basedir/results/$1/armadillo.log
  exit 1
fi

# Zip and copy the results to s3
rm -f $basedir/results/zip/$1_result.zip
rm -f $basedir/results/zip/$1_report.zip
zip -rj $basedir/results/zip/$1_report $basedir/results/$1/*.tsv $basedir/results/$1/*.pdf $basedir/results/$1/*.xlsx $basedir/results/$1/*.log $basedir/results/$1/phx_output/multiqc/multiqc_report.html $basedir/results/$1/phx_output/*.tsv $basedir/results/$1/phx_output/*.xlsx
aws s3 cp $basedir/results/zip/$1_report.zip $aws_bucket/ARLN/REPORT/$1_report.zip --region us-gov-west-1

zip -r $basedir/results/zip/$1_result $basedir/results/$1
aws s3 cp $basedir/results/zip/$1_result.zip $aws_bucket/ARLN/ANALYSIS_RESULTS/$1_result.zip --region us-gov-west-1

rm $basedir/results/zip/$1_*.zip 
rm $basedir/reads/zip/$1.zip 
