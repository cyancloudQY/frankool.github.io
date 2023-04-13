# Detailed Notes about Including Repetive Elements in RNA-seq to counts workflow (TEtranscripts)
1. Step-by-Step illustration of the upstream analysis of RNA-seq data, while including the quantification of the expression of transposable elements. 
2. Didn't include quality control in the workflow, need to add it afterwards. 
3. Remember, should never remove duplicates in RNA-seq data. 


## 1. Download Required Software
### TEtranscripts

```bash
#Install anaconda
wget https://mirrors.bfsu.edu.cn/anaconda/archive/Anaconda3-2022.10-Linux-x86_64.sh --no-check-certificate
bash Anaconda3-2022.10-Linux-x86_64.sh 

# Install pysam
pip3 install pysam

# This will install the package in the user directory, need to add the path into the environment
git clone https://github.com/mhammell-laboratory/TEtranscripts
cd TEtranscripts
python3 setup.py install --user
export PATH="$HOME/.local/bin:$PATH"

# Get the curated gtf file for TE
wget https://labshare.cshl.edu/shares/mhammelllab/www-data/TEtranscripts/TE_GTF/GRCh38_GENCODE_rmsk_TE.gtf.gz


# Get the Genecode GTF file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.annotation.gtf.gz

```

### STAR
```bash
wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
tar -xzf 2.7.10b.tar.gz
cd STAR/source
make STAR

# Human genome fasta file
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.p13.genome.fa.gz

# Generate genome index directory. (TAKES HOURS)
/home/yojeep/Desktop/Bioinfo/app/STAR/STAR-2.7.10b/bin/Linux_x86_64/STAR  \
--runThreadN 6 \
--runMode genomeGenerate  \
--genomeDir ./  \
--genomeFastaFiles /home/yojeep/Desktop/Bioinfo/reference/Genecode/GRCh38.p13.genome.fa  \
--sjdbGTFfile /home/yojeep/Desktop/Bioinfo/reference/Genecode/gencode.v43.annotation.gtf 
```


### Aspera
A super fast tool to download sequencing data. 
```bash
wget -c https://ak-delivery04-mul.dhe.ibm.com/sar/CMA/OSA/09ff1/0/ibm-aspera-connect-3.11.1.58-linux-g2.12-64.tar.gz

tar xzvf ibm-aspera-connect-3.11.1.58-linux-g2.12-64.tar.gz
bash ibm-aspera-connect-3.11.1.58-linux-g2.12-64.sh

bash ibm-aspera-connect-3.11.1.58-linux-g2.12-64.sh

vi ~/.bashrc
export PATH="$PATH:$HOME/.aspera/connect/bin"

```


## 2. Download + Align + Count
Download data table (tsv) including the asepra location from EMBL-EBI using the BioProject ID or SRA ID

```python
import pandas as pd

## ebi
ebi_file = pd.read_csv('./filereport_read_run_PRJNA505280_tsv.txt', sep = '\t')

fastq_link = ebi_file['fastq_aspera']
fastq_link = fastq_link.to_list()

sep_link = []
for i in fastq_link:
    a,b = i.split(';')
    sep_link.append(a)
    sep_link.append(b)

aspera_l = 'ascp -k 2 -QT -l 20m -P33001 -i $HOME/.aspera/connect/etc/asperaweb_id_dsa.openssh era-fasp@'
f = open('./ascp_list.sh','w')

for i in sep_link:
    command = aspera_l + i + ' ./ \n'
    f.writelines(command)

f.close()
print(sep_link)
```


```bash
# Download
bash ascp_list.sh
```

Align the fasta using STAR, need to set the parameter to include the multi-mapped reads.
```bash
/home/yojeep/Desktop/Bioinfo/app/STAR/STAR-2.7.10b/bin/Linux_x86_64/STAR \
--genomeDir /home/yojeep/Desktop/Bioinfo/reference/Genecode/STAR \
--readFilesIn /home/yojeep/Desktop/Bioinfo/workplace/PRJNA505280/fasta/SRR8181359_1.fastq.gz /home/yojeep/Desktop/Bioinfo/workplace/PRJNA505280/fasta/SRR8181359_2.fastq.gz \
--outFileNamePrefix SRR8181359 \
--sjdbScore 2 --runThreadN 6 \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate \
--winAnchorMultimapNmax 100 \
--outFilterMultimapNmax 100
```

OK, good to count the bam file using TEcount. Remember to add proper prefix in order for smooth downstream analysis.
```
for i in *.bam
do
TEcount                      \
	-b $i            \
	--TE /home/yojeep/Desktop/Bioinfo/reference/TEtranscripts/GRCh38_GENCODE_rmsk_TE.gtf         \
	--GTF /home/yojeep/Desktop/Bioinfo/reference/Genecode/gencode.v43.annotation.gtf \
	--sortByPos \
	--project "${i%%.*}"
done
```

**Some Thoughts**: Should write above analysis including quality control into a scripts that process one sample at a time. Need to consider the file location.


## 3. Aggregate raw counts using R
Extract GTF info for gene-to-symbol conversion.
```bash

gtf="gtf_file_name.gtf"
### gene_id to gene_name
grep 'gene_id' $gtf | awk -F 'gene_id \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_id_tmp
grep 'gene_id' $gtf | awk -F 'gene_name \"' '{print $2}' |awk -F '\"' '{print $1}' >gene_name_tmp
paste gene_id_tmp gene_name_tmp >last_tmp
uniq last_tmp >g2s_gencode.txt
rm *_tmp

```




Use R to aggragate the raw counts.

```R
library(tidyverse)
library(DESeq2)


## All the raw count file should be located in a folder with proper 
## name that specify their origins of samples

## Merge the raw counts
files <- list.files(path = './TE_count/')

merge.data <- read.table(paste('./TE_count/', files[1], sep = ''), row.names = 1, header = TRUE)

for (i in 2:length(files)){
  new.data <- read.table(paste('./TE_count/', files[i], sep = ''), row.names = 1, header = TRUE)
  merge.data <- cbind(merge.data, new.data)
}

## Extract (the variable name is wrong, should be all transposons)
all_gene <- rownames(merge.data)
LINE <- all_gene[grep(':', all_gene)]
LINE_g2s <- data.frame(ensemble = LINE, symbol = LINE)

## Gene symbol
g2s <- read.table('./g2s_human_gencode.txt')
colnames(g2s) <- c('ensemble', 'symbol')
g2s <- rbind(g2s, LINE_g2s)


symbol <- g2s[match(rownames(merge.data),g2s$ensemble),"symbol"]
table(duplicated(symbol))

counts <- aggregate(merge.data, by=list(symbol), FUN=sum)
counts <- column_to_rownames(counts,'Group.1')


## DEG
group_info <- factor(c(rep('SLE', 3), rep('control', 3)), levels = c('SLE', 'control'))

DESeq_object <- DESeqDataSetFromMatrix(counts, colData = DataFrame(group_info), design= ~group_info)
DESeq_object <- DESeq(DESeq_object)
DESeq_results <- as.data.frame(results(DESeq_object))

write.table(DESeq_results, 'a.csv', sep = '\t')


DESeq_results_TE <- filter(DESeq_results, rownames(DESeq_results) %in% LINE_g2s$symbol)

DESeq_results_gene <- filter(DESeq_results, !rownames(DESeq_results) %in% LINE_g2s$symbol)

```