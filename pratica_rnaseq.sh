
mkdir -p PRATICA_RNASEQ
cd PRATICA_RNASEQ
scp -P 1000 jorge.munoz@bioinfo.cena.usp.br:/Storage/data1/jorge.munoz/PRACTICA_RNASEQ/data/tiny/data.zip .
mkdir data
cd data
unzip ./../data.zip
cd ..
rm data.zip





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

# create tx2gene
gunzip -c ./references/gencode.v38.transcripts.fa.gz \
| grep ">"| cut -f 1,2 -d '|' \
|sed 's/|/\t/'|sed 's/>//' > tx2gene.txt
# run fastqc
mkdir -p QC
conda activate fastqc_env
for file in $(ls -1 data/*.fastq); do fastqc --noextract \
--threads 1 --nogroup -o QC $file; done
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
