https://linuxize.com/post/how-to-install-r-on-ubuntu-20-04/

#mkdir -p PRATICA_RNASEQ
#cd PRATICA_RNASEQ
cd 
scp -P 1000 jorge.munoz@bioinfo.cena.usp.br:/Storage/data1/jorge.munoz/PRACTICA_RNASEQ/PRATICA_RNASEQ.zip .
#scp -P 1000 jorge.munoz@bioinfo.cena.usp.br:/Storage/data1/jorge.munoz/PRACTICA_RNASEQ/data/tiny/data.zip .
#mkdir data
#cd data
unzip PRATICA_RNASEQ.zip
#cd ..
rm PRATICA_RNASEQ.zip

apt update -y && apt upgrade -y
# Install r
# update indices
apt update -qq
# install two helper packages we need
apt install --no-install-recommends software-properties-common dirmngr
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc 
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
apt install --no-install-recommends r-base
# Download Rstudio
wget https://download1.rstudio.org/desktop/bionic/amd64/rstudio-2022.02.3-492-amd64.deb
# install Rstudio
dpkg -i rstudio-2022.02.3-492-amd64.deb


sudo apt install libxml2-dev libcurl4-openssl-dev libssl-dev


#######
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("tximeta")
install.packages(c("pheatmap","viridis","ggplot2"))


# go wd
cd ~/PRACTICA_RNASEQ
mkdir -p prefetch
# Download accessions
for i in $(cut -f1 -d"," metadata.csv | grep SRR)
do
	prefetch -v $i -O ./prefetch
done

conda create -n parallel_fastq_dump
conda activate parallel_fastq_dump	
conda install -c bioconda parallel-fastq-dump 

for i in $(cut -f1 -d"," metadata.csv | grep SRR)
do
	parallel-fastq-dump --threads 4 --outdir ./data/ --split-files --gzip --sra-id ./prefetch/${i}.sra
done
# download references
cd PRATICA_RNASEQ
mkdir references
cd references
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz
# install salmon
wget https://github.com/COMBINE-lab/salmon/releases/download/v1.8.0/salmon-1.8.0_linux_x86_64.tar.gz
tar -xvzf salmon-1.8.0_linux_x86_64.tar.gz
rm salmon-1.8.0_linux_x86_64.tar.gz
sudo echo export PATH=$PATH:~/salmon-1.8.0_linux_x86_64/bin >> ~/.bashrc
source ~/.bashrc
# create salmon index
cd ~/PRATICA_RNASEQ
salmon index -t ./references/gencode.v38.transcripts.fa.gz -i salmon_index --threads 2


# run fastqc1
mkdir -p QC1
conda activate fastqc_env
for file in $(ls data/*.fastq); do fastqc --noextract \
--threads 6 --nogroup -o QC1 $file; done
conda deactivate


# run bbduk
mkdir -p bbduk
conda activate bbduk_env
for sample in $(cut -f1 -d"," metadata.csv | grep SRR)
do
file_R1=${sample}_1.fastq 
file_R2=${sample}_2.fastq
bbduk.sh in1=./data/${file_R1} in2=./data/${file_R2} out1=./bbduk/${file_R1} out2=./bbduk/${file_R2} minlength=50 qtrim=w trimq=20
done
conda deactivate

# run fastqc2
mkdir -p QC2
conda activate fastqc_env
for file in $(ls data/*.fastq); do fastqc --noextract \
--threads 6 --nogroup -o QC2 $file; done
conda deactivate

# run salmon
mkdir -p Quantification
cd Quantification
for sample in $(cut -f1 -d"," ./../metadata.csv | grep SRR); do \
file_R1=${sample}_1.fastq; \
file_R2=${sample}_2.fastq; \
salmon quant --libType A --threads 2 --index \
./../salmon_index/ \
--validateMappings --seqBias --posBias --softclip \
-1 ./../bbduk/${file_R1} \
-2 ./../bbduk/${file_R2} \
-o ${sample};
done


## BLAST PARA FILTAR DATOS DEL CROMOSOMA 1
cd data
mkdir blast
export BLASTDB=/home/j/PRACTICA_RNASEQ/data/blast
cd blast
makeblastdb -in ./../../salmon/chr1.fna -dbtype 'nucl' -out myDB
cd ..
blastn -db myDB -query $(zcat ./../data/SRR1039508_1.fastq.gz)
