#!/bin/bash
# Define the list of files to copy from S3
FILES_TO_COPY=(
	"call_run_pipeline_prod.sh"
    "call_run_pipeline_test.sh"
	"run_pipeline_prod.sh"
    "run_pipeline_test.sh"
)

# S3 bucket path
S3_BUCKET="s3://430118851772-bioinformatics-code/scripts"

# Copy each file from S3
for file in "${FILES_TO_COPY[@]}"; do
	aws s3 cp "$S3_BUCKET/$file" .
    # Set permissions for copied scripts
    chmod 755 "$file"
done

