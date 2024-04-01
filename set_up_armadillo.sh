git clone --branch bfx https://github.com/TX-DSHS/Armadillo_pipeline.git
mkdir reads
mkdir reads/zip
mkdir scripts
mkdir results
mkdir results/zip
mkdir cluster
aws s3 cp --recursive s3://804609861260-bioinformatics-infectious-disease/pipeline/kraken2_db/k2_standard_08gb_20230605 ./ --region="us-gov-west-1"
aws s3 cp --recursive s3://804609861260-bioinformatics-infectious-disease/pipeline/ARLN/Phoenix/miniconda3 ./ --region="us-gov-west-1"
aws s3 cp --recursive s3://804609861260-bioinformatics-infectious-disease/pipeline/ARLN/Phoenix/phoenix ./ --region="us-gov-west-1"
# modify aws_bucket.txt file
