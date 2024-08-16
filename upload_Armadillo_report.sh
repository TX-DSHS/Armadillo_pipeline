basedir=$PWD
rm -r $basedir/results/$1/work
rm $basedir/results/zip/$1_report.zip
zip -rj $basedir/results/zip/$1_report $basedir/results/$1/*.tsv $basedir/results/$1/*.pdf $basedir/results/$1/*.xlsx $basedir/results/$1/*.log $basedir/results/$1/phx_output/multiqc/multiqc_report.html $basedir/results/$1/phx_output/*.tsv $basedir/results/$1/phx_output/*.xlsx
aws s3 cp $basedir/results/zip/$1_report.zip s3://804609861260-bioinformatics-infectious-disease/ARLN/REPORT/$1_report.zip --region us-gov-west-1
zip -r $basedir/results/zip/$1_result $basedir/results/$1
aws s3 cp $basedir/results/zip/$1_result.zip s3://804609861260-bioinformatics-infectious-disease/ARLN/ANALYSIS_RESULTS/$1_result.zip --region us-gov-west-1
