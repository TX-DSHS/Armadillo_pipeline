# Armadillo_pipeline
Texas DSHS Armadillo pipeline is designed to process ARLN samples using CDC PhoeNix pipeline:
https://github.com/CDCgov/phoenix

## The pipeline can be installed in /bioinformatics/ partition of AWS EC2 by:
```bash
git clone https://github.com/TX-DSHS/Armadillo_pipeline.git -b bfx

# Sometimes if there is "\r" or EOF error
git config core.autocrlf false
git rm --cached -r .         # Don’t forget the dot at the end
git reset --hard
```
## Install conda
```bash
curl -sL \
  "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh" > \
  "Miniconda3.sh"
bash Miniconda3.sh
```
Do you accept the license terms? [yes|no]
>>> yes

Miniconda3 will now be installed into this location:
/home/dnalab/miniconda3

  - Press ENTER to confirm the location
  - Press CTRL-C to abort the installation
  - Or specify a different location below

[/home/dnalab/miniconda3] >>> /bioinformatics/Armadillo_pipeline/miniconda3

```bash
rm Miniconda3.sh
```

## Create a conda environment installing Singularity and nextflow:
```bash
source /bioinformatics/Armadillo_pipeline/miniconda3/etc/profile.d/conda.sh
conda create -n nextflow -c conda-forge -c bioconda openjdk==11.0.20 singularity nextflow
conda install -n nextflow --file python_requirements.txt
```

## Download Kraken2 databases and phoenix scripts
```bash
aws s3 cp --recursive s3://804609861260-bioinformatics-infectious-disease/pipeline/kraken2_db/k2_standard_08gb_20230605 ./ --region="us-gov-west-1"
aws s3 cp --recursive s3://804609861260-bioinformatics-infectious-disease/pipeline/ARLN/Phoenix/phoenix ./ --region="us-gov-west-1"
sudo chmod +x /bioinformatics/Armadillo_pipeline/phoenix/bin/*
```
## Contributing

Pull requests are welcome. For major changes, please open an issue first
to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License

[MIT](https://choosealicense.com/licenses/mit/)
