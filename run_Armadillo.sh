basedir="/home/dnalab"
aws s3 cp s3://804609861260-bioinformatics-infectious-disease/ARLN/RAW_RUNS/$1.zip $basedir/reads/zip --region us-gov-west-1
mkdir $basedir/reads/$1
unzip -j $basedir/reads/zip/$1.zip -d $basedir/reads/$1
mkdir $basedir/results/$1
cd $basedir/results/$1

nextflow run /pipeline/ARLN/Grandeur/grandeur.nf -profile singularity -c /pipeline/ARLN/Grandeur/configs/AMR.config --reads $basedir/reads/$1 --outdir $basedir/results/$1

rm -r $basedir/results/$1/work
rm -r $basedir/results/$1/shuffled

python3 /home/jiel/bin/armadillo/armadillo.py grandeur_results.tsv

rm -f $basedir/results/zip/$1_result.zip
rm -f $basedir/results/zip/$1_report.zip
zip -rj $basedir/results/zip/$1_report $basedir/results/$1/*.tsv $basedir/results/$1/*.pdf $basedir/results/$1/*.xlsx $basedir/results/$1/*.log
zip -r $basedir/results/zip/$1_result $basedir/results/$1
aws s3 cp $basedir/results/zip/$1_report.zip s3://804609861260-bioinformatics-infectious-disease/ARLN/REPORT/$1_report.zip --region us-gov-west-1
aws s3 cp $basedir/results/zip/$1_result.zip s3://804609861260-bioinformatics-infectious-disease/ARLN/ANALYSIS_RESULTS/$1_result.zip --region us-gov-west-1
