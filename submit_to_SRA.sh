mkdir /home/dnalab/reads/$1/sra
cp -d /home/dnalab/reads/$1/20*.fastq.gz /home/dnalab/reads/$1/sra/
/home/jiel/.aspera/cli/bin/ascp -i /home/dnalab/sra/aspera.openssh -QT -l100m -k1 -d /home/dnalab/reads/$1/sra subasp@upload.ncbi.nlm.nih.gov:uploads/lab.microbiology_dshs.state.tx.us_rJiZeQDA/$1
rm -r /home/dnalab/reads/$1/sra
