# Run Armadillo pipeline
# Author: Jie.Lu@dshs.texas.gov
version="2.0-10/2/2023"

basedir="/home/dnalab"
rm -rf $basedir/results/$1
mkdir -p $basedir/results/$1
echo "Running Armadillo version: "$version > $basedir/results/$1/armadillo.log

# # Copy and unzip the fastq files from s3
aws s3 cp s3://804609861260-bioinformatics-infectious-disease/ARLN/RAW_RUNS/$1.zip $basedir/reads/zip --region us-gov-west-1
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

# Run PhoeNix pipeline
reads_path=$basedir/reads/$1
out_path=$basedir/results/$1
$basedir/phoenix/bin/create_samplesheet.sh $reads_path > $reads_path/samplesheet.csv
source $basedir/miniconda3/etc/profile.d/conda.sh
conda activate nextflow
cd $out_path
export NXF_SINGULARITY_CACHEDIR=/singularity/Phoenix
nextflow run cdcgov/phoenix -r v2.0.2 -profile singularity -entry PHOENIX --input $reads_path/samplesheet.csv --kraken2db $basedir/kraken2_db/k2_standard_08gb_20230605 --output $out_path
conda deactivate
rm -r $out_path/work
rm $basedir/results/zip/$1.zip

# Run Armadillo pipeline
python3 /home/jiel/bin/armadillo/armadillo_phoenix.py -r $1 >> $basedir/results/$1/armadillo.log

# Zip and copy the results to s3
rm -f $basedir/results/zip/$1_result.zip
rm -f $basedir/results/zip/$1_report.zip
#zip -rj $basedir/results/zip/$1_report $basedir/results/$1/*.tsv $basedir/results/$1/*.pdf $basedir/results/$1/*.xlsx $basedir/results/$1/*.log $basedir/results/$1/results/multiqc/multiqc_report.html $basedir/results/$1/results/Phoenix_Summary.tsv
#zip -rj $basedir/results/zip/$1_report $basedir/results/$1/*.log $basedir/results/$1/results/multiqc/multiqc_report.html $basedir/results/$1/results/Phoenix_Summary.tsv

#aws s3 cp $basedir/results/zip/$1_report.zip s3://804609861260-bioinformatics-infectious-disease/ARLN/REPORT/$1_report.zip --region us-gov-west-1
#zip -r $basedir/results/zip/$1_result $basedir/results/$1
#aws s3 cp $basedir/results/zip/$1_result.zip s3://804609861260-bioinformatics-infectious-disease/ARLN/ANALYSIS_RESULTS/$1_result.zip --region us-gov-west-1

rm $basedir/results/zip/$1_*.zip 
rm $basedir/reads/zip/$1.zip 
