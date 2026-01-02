#!/bin/bash
# This script calls run_pipeline.sh with 3 parameters

# Example parameters (replace with actual values as needed)
RUN_ID="AR_251124C_VH00729"
S3_BUCKETNAME="430118851772-bioinformatics-infectious-disease"
SNS_TOPIC_ARN="arn:aws:sns:us-east-1:430118851772:dshs-bioinfo-bacteria-notification"

# Call the pipeline script
"$(dirname "$0")/run_pipeline_prod.sh" "$RUN_ID" "$S3_BUCKETNAME" "$SNS_TOPIC_ARN"
