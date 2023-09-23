# Hadins
Hadins: Insertion detection based on deep learning using high accuracy data

## Installation
### Requirements
* python 3.9, numpy, pandas, Matplotlib, TensorFlow 2.7, pysam
### 1. Create a virtual environment  
```
#create
conda create -n Hadins python=3.9
#activate
conda activate Hadins
#deactivate
conda deactivate
```   
### 2. clone LSnet
* After creating and activating the Hadins virtual environment, download Hadins from github:
```　 
https://github.com/cyq1314woaini/Hadins.git
cd Hadins
```
### 3. Install 
```　
conda activate Hadins
conda install numpy, pandas, Matplotlib, TensorFlow 2.7, pysam
```
## Usage
### 1.Generate Features
```　 
python Hadins.py generate_feature lbamfile sbamfile output process includecontig\n
lbamfile is the path of the alignment file about the reference and the long read set\n
sbamfile is the path of the alignment file about the reference and the short read set\n
output is a folder which is used to store features\n
process is the number of processes to use\n
includecontig is the list of contig to preform detection.(default: [], all contig are used) \n
``` 
### 2.Call Insertion 
```　 
python Hadins.py weight datapath lbamfile outvcfpath support includecontig\n
weight is the path of the model weight\n
datapath is a folder which is used to store features\n
lbamfile is the path os the alignment file about the reference and the long read set\n
outvcfpath is the path of vcf file\n
support is min support reads\n
includecontig is the list of contig to preform detection.(default: [], all contig are used)\n 
```  
## Tested data 
The data can be downloaded from:  

### HG002 long read data
```
https://ftp.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_MtSinai_NIST/Baylor_NGMLR_bam_GRCh37/HG002_PB_70x_RG_HP10XtrioRTG.bam
```
### HG002 short read data
```
https://ftp.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.hs37d5.2x250.bam
```   
### HG002 ccs data
``` 
https://ftp.ncbi.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/alignment/HG002.Sequel.15kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
 ```  
### NA19240 long read data
```
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20160905_smithm_pacbio_aligns/NA19240_bwamem_GRCh38DH_YRI_20160905_pacbio.bam
```   
### NA19240 short read data
```
http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/YRI/NA19240/high_cov_alignment/NA19240.alt_bwamem_GRCh38DH.20150715.YRI.high_coverage.cram
```   
