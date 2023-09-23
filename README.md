# Hadins
Hadins: Insertion detection based on deep learning using high accuracy data


## Installation
### Requirements
* python 3.9, numpy, pandas, Matplotlib, TensorFlow 2.7, pysam

### 1. Create a virtual environment  
'''
#create
conda create -n Hadins python=3.9
#activate
conda activate Hadins
'''

### 2. clone Hadins
'''
cd https://github.com/cyq1314woaini/Hadins.git
cd Hadins
'''

### 3.Install
'''
conda activate LSnet
conda install numpy, pandas, Matplotlib, TensorFlow 2.7, pysam
'''

## Usage
### 1.Generate Feature
'''
python Hadins.py generate_feature lbamfile sbamfile output process includecontig
lbamfile is the path of the alignment file about the reference and the long read set
sbamfile is the path of the alignment file about the reference and the short read set
output is a folder which is used to store features
process is the number of processes to use
includecontig is the list of contig to preform detection.(default: [], all contig are used)
'''

### 2.Call Insertion
'''
python Hadins.py weight datapath lbamfile outvcfpath support includecontig\n
weight is the path of the model weight
datapath is a folder which is used to store features
lbamfile is the path os the alignment file about the reference and the long read set
outvcfpath is the path of vcf file
support is min support reads
includecontig is the list of contig to preform detection.(default: [], all contig are used)
'''
