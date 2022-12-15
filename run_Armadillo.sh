# Run Armadillo pipeline
# Author: Jie.Lu@dshs.texas.gov
version="1.0-10/11/2022"

basedir="/home/dnalab"
mkdir -p $basedir/results/$1
echo "Running Armadillo version %s" $version > $basedir/results/$1/armadillo.log

# Copy and unzip the fastq files from s3
aws s3 cp s3://804609861260-bioinformatics-infectious-disease/ARLN/RAW_RUNS/$1.zip $basedir/reads/zip --region us-gov-west-1
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

# Run Grandeur pipeline
cd $basedir/results/$1
nextflow run /pipeline/ARLN/Grandeur/grandeur.nf -profile singularity -c /pipeline/ARLN/Grandeur/configs/AMR.config --reads $basedir/reads/$1 --outdir $basedir/results/$1
rm -r $basedir/results/$1/work
rm -r $basedir/results/$1/shuffled
rm $basedir/results/zip/$1.zip
# Run Armadillo pipeline
python3 /home/jiel/bin/armadillo/armadillo.py -r $1

# Zip and copy the results to s3
rm -f $basedir/results/zip/$1_result.zip
rm -f $basedir/results/zip/$1_report.zip
zip -rj $basedir/results/zip/$1_report $basedir/results/$1/*.tsv $basedir/results/$1/*.pdf $basedir/results/$1/*.xlsx $basedir/results/$1/*.log
aws s3 cp $basedir/results/zip/$1_report.zip s3://804609861260-bioinformatics-infectious-disease/ARLN/REPORT/$1_report.zip --region us-gov-west-1
zip -r $basedir/results/zip/$1_result $basedir/results/$1
aws s3 cp $basedir/results/zip/$1_result.zip s3://804609861260-bioinformatics-infectious-disease/ARLN/ANALYSIS_RESULTS/$1_result.zip --region us-gov-west-1

rm $basedir/results/zip/$1_*.zip 
rm $basedir/reads/zip/$1.zip 

