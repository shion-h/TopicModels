# TopicModels
# Abstract
This C++ code implements Bayesian inference algorithm of topic models.
# Installation
    cmake . -DCMAKE_CXX_COMPILER=g++
    make
# Requirements
The boost library is required.  

Regarding python scripts, you should install numpy, scipy, and pandas.  
# Usage
## Latent Dirichlet allocation
 LDA [BOW file] [-options] 

Options:  
  -h [ --help ]              show help  
  -k [ --nmtp ] arg (=10)    the number of topics  
  -d [ --cdrt ] arg (=0.001) convergence ditermination rate  
  -o [ --otpt ] arg (=./)    directory name for output  
  
## Author topic model
 ATM [BOW file] [author file] [-options] 

Options:  
  -h [ --help ]              show help  
  -k [ --nmtp ] arg (=10)    the number of topics  
  -d [ --cdrt ] arg (=0.001) convergence ditermination rate  
  -o [ --otpt ] arg (=./)    directory name for output  

## Supervised topic model
 sLDA [BOW file] [label file] [-options] 

Options:  
  -h [ --help ]              show help  
  -k [ --nmtp ] arg (=10)    the number of topics  
  -d [ --cdrt ] arg (=0.001) convergence ditermination rate  
  -o [ --otpt ] arg (=./)    directory name for output  

## Hierarchical Dirichlet process
 HDP [BOW file] [-options] 

Options:  
  -h [ --help ]             show help  
  -o [ --otpt ] arg (=./)   directory name for output  
  -n [ --itnm ] arg (=1000) the number of iteration  
  -i [ --intr ] arg (=5)    sampling interval  
  -b [ --bnin ] arg (=500)  burn-in term  
  -a [ --alph ] arg (=100)  alpha value  
  -c [ --beta ] arg (=100)  beta value  
  -g [ --gamm ] arg (=100)  gamma value  
# Input
  An input file is required to be frequency matrix format. See data/test.csv.  
# Output
## LDA
  1)topic distribution(theta.csv)  
  2)word distribution(phi.csv)  
  3)hyper parameters(alpha.csv, beta.csv)  
  4)word list(wordList)  
  ( 5)variational lower bound(variationalLowerBound))  
# Scripts
  You can also use python scripts. These scripts calculate some values for each the number of topics.  
## Usage
### Perplexity
calcldappl.py [-h] [-s KSTEP] [-r RATE] [-i ITERATION]  
                     BOW_filename k_min k_max  
  
positional arguments:  
  - BOW_filename:          Bag Of Words filename  
  - k_min:                 minimum number of topics for calculation of perplexity  
  - k_max:                 max number of topics for calculation of perplexity  
  
optional arguments:  
  - -h, --help:            show this help message and exit  
  - -s KSTEP, --kstep KSTEP: step number of topics for calculation of perplexity  
  - -r RATE, --rate RATE:  rate of train data  
  - -d CONV_DET, --conv_det CONV_DET: convergence determination
  
### WBIC
calcldawbic.py [-h] [-s KSTEP] [-i ITERATION] BOW_filename k_min k_max  
  
positional arguments:  
  - BOW_filename:          Bag Of Words filename  
  - k_min:                 minimum number of topics for calculation of wbic  
  - k_max:                 max number of topics for calculation of wbic  
  
optional arguments:  
  - -h, --help:            show this help message and exit  
  - -s KSTEP, --kstep KSTEP: step number of topics for calculation of wbic  
  - -d CONV_DET, --conv_det CONV_DET: convergence determination
  
### Maximized variational lower bound
deriveldavlb.py [-h] [-s KSTEP] BOW_filename k_min k_max  
  
positional arguments:  
  - BOW_filename:          Bag Of Words filename  
  - k_min:                 minimum number of topics for calculation of maximized variational lower bound  
  - k_max:                 max number of topics for calculation of maximized variational lower bound  
  
optional arguments:  
  - -h, --help:            show this help message and exit  
  - -s KSTEP, --kstep KSTEP: step number of topics for calculation of maximized variational lower bound  
  - -d CONV_DET, --conv_det CONV_DET: convergence determination
  - -o OUTPUT_DIR, --output_dir OUTPUT_DIR: output directory
## Output
These scripts output results to /output/(time stamp) directory by default. (time stamp) indecates the date information in MMDD_DAY_hh_mm_ss format.  
# License
This software is released under the MIT License, see LICENSE.txt.  
# How to run for modeling taxonomic profiles [Hosoda et al., Microbiome, 2020](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00864-3)
1. Requirements
    1. Taxonomic profile dataset (Same form as data/test.csv)
        - This is one CSV format file (seperated by commas).
        - This has D+1 rows and V+1 columns, where D is the number of samples and V is the number of different microbes.
        - The 1st row indicates the names of microbes.
        - The 1st column indicates the names of samples.
        - Each element is a disrete number.
2. Work on a terminal
    ```
    cd /path/to/TopicModels/root
    ./LDA /path/to/your/taxonomic/profile/dataset.csv -k 4 -d 1.0e-6 -o /path/to/your/output/directory
    ```
    - The number of assemblages (four in the example above) can be changed as needed.
3. Work in the downstream
    1. You can see theta.csv and phi.csv in /path/to/your/output/directory
    2. theta.csv is the assemblage distributions for each samples.
        - This has D rows and K columns, where D is the number of samples and K is the number of assemblages.
    3. phi.csv is the microbe distributions for each assemblage.
        - This has K rows and V columns, where K is the number of assemblages and V is the number of different microbes.
