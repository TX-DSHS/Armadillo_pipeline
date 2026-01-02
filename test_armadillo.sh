source /bioinformatics/Armadillo_pipeline/miniconda3/etc/profile.d/conda.sh
conda activate nextflow
python3 ../../armadillo_phoenix.py -r AR_251219_M04955 -d /bioinformatics/Armadillo_pipeline -a s3://430118851772-bioinformatics-infectious-disease -c /ReferenceGeneCatalog_3.12_20240205.txt
