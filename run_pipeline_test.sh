#!/bin/bash

# -----------------------------------------------------------------------------
# Script: run_pipeline_test.sh
# Description:
#   This script tests the setup for the Armadillo pipeline run in a test environment. 
#   It uses the root folder /bioinformatics/Armadillo/test to create all folders.
#   
#   It sets up the required directory structure, downloads input data from S3,
#   manages logs, and sends notifications via AWS SNS. It does not run the full pipeline but is intended for testing the setup.
#
# Required Input Parameters:
#   1. run_id         - Unique identifier for the pipeline run.  This is the filename prefix. (e.g. AR_251124C_VH00729)
#   2. s3_bucketname  - Name of the AWS S3 bucket containing input data (e.g. 430118851772-bioinformatics-infectious-disease)
#   3. sns_topicarn   - ARN of the AWS SNS topic for notifications (e.g. arn:aws:sns:us-east-1:430118851772:dshs-bioinfo-bacteria-notification)
#
# Usage:
#   bash run_pipeline_test.sh <run_id> <s3_bucketname> <sns_topicarn>
#
# Example:
#   bash run_pipeline_test.sh AR_251124C_VH00729 430118851772-bioinformatics-infectious-disease arn:aws:sns:us-east-1:430118851772:dshs-bioinfo-bacteria-notification
# -----------------------------------------------------------------------------



# Global variable definitions
version="v1.0.0"
basedir=$PWD/test # for testing
# first log file to track all the setting before the run begins
logfile="$basedir/run_pipeline.log"

# if the logfile already exists, remove it
[ -f "$logfile" ] && rm "$logfile"
####################################################################
# Function definitions
####################################################################
# Function to send a message to SNS
send_sns_message() {
  local message="$1"
  aws sns publish --topic-arn "$sns_topicarn" --message "$message"
  if [ $? -eq 0 ]; then
    echo "SNS message sent: $message" >> "$logfile"
  else
    echo "Failed to send SNS message: $message" >> "$logfile"
  fi
}


# Function to handle failures by storing it in the log file and sending SNS message and zipping log files and storing them in S3
handle_failure() {
  local message="$1"
  # Store in log file
  echo "$message" >> "$logfile"
  # Send the SNS message
  send_sns_message "$message"
  # zip the log file
  zip -rj $basedir/results/zip/${1}_report $basedir/results/$1/*.log
  # Copy the log file to S3
  #aws s3 cp $basedir/results/zip/${1}_report.zip $aws_bucket/ARLN/REPORT/${1}_report.zip
  # Exit the script
  exit 1
}

# Function to validate run_id length.
# This is an example of how the validation can be done in the main script.
validate() {
  local run_id="$1"
  # if [ ${#run_id} -le 12 ]; then
  #   message ="Error: run_id '$run_id' must be more than 12 characters."
  #   send_sns_message "$message"
  #   echo "$message" >> "$logfile"
  #   exit 1
  # fi
}


# Ensure the base directory is created
mkdir -p $basedir

echo "Script Arguments:" >> "$logfile"
for arg in "$@"; do
  echo "Arg: $arg" >> "$logfile"
done

if [ $# -ne 3 ]; then
  echo "Usage: bash run_pipeline.sh run_id s3_bucketname sns_topicarn" >> "$logfile"
  echo "Exiting" >> "$logfile"
  exit 1
fi


date >> $logfile


# runid  is the filename prefix (example AR_251124C_VH00729)
run_id=$1
# S3 bucketname (e.g. 430118851772-bioinformatics-infectious-disease)
s3_bucketname=$2
aws_bucket="s3://$s3_bucketname"
# SNS topic ARN to send notification (e.g. arn:aws:sns:us-east-1:430118851772:dshs-bioinfo-bacteria-notification) 
sns_topicarn=$3


# example of how to use validatition on the run_id prior to running the pipeline
validate "$run_id"


# make necessary root directories if they don't already exist
mkdir -p $basedir/results
mkdir -p $basedir/reads
mkdir -p $basedir/results/zip
mkdir -p $basedir/reads/zip
mkdir -p $basedir/cluster


# clean up any previous run data
rm -rf $basedir/results/$1
mkdir -p $basedir/results/$1

rm -f $basedir/results/zip/$1_result.zip
rm -f $basedir/results/zip/$1_report.zip


# create the log file for this run
run_logfile="$basedir/results/$1/armadillo.log"
echo "Run Logfile: $run_logfile" >> "$logfile"


####################################################################
# Pipeline Starting Point
####################################################################
echo "Pipeline Version: "$version >> $run_logfile
date > $run_logfile
echo "Run ID: $1" >> $run_logfile
echo "AWS S3 Bucket: $aws_bucket" >> $run_logfile
echo "SNS Topic ARN: $sns_topicarn" >> $run_logfile
# Send SNS notification that the pipeline run has started
message="Pipeline Started Successfully on $(date) with run_id: $1"
echo "$message" >> $run_logfile
send_sns_message "$message"



# Copy and unzip the fastq files from s3
s3_object="$aws_bucket/ARLN/RAW_RUNS/$1.zip"
aws s3 cp $s3_object $basedir/reads/zip
# if upload failed, upload the log file to aws, exit
if [ $? -ne 0 ]; then
  message="Failed to download the fastq files from s3_object $s3_object for run_id: $1"
  handle_failure "$message"
fi


# Stop Testing here.




message="Pipeline Finished Successfully on $(date) with run_id: $1"
echo "$message" >> $run_logfile
send_sns_message "$message"



