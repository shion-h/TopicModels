# TopicModels
# Abstract
These tools were developed to estimate parameters (including hyper parameters) of some topic models.  
# Installation
    cmake . -DCMAKE_CXX_COMPILER=g++
    make
# Requirements
Installing the boost library is required for using this tools.  

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
  This program outputs the following:  
## LDA
  1)topic distribution(theta.csv)  
  2)word distribution(phi.csv)  
  3)hyper parameters(alpha.csv, beta.csv)  
  4)word list(wordList)  
  ( 5)vatiational lower bound(variationalLowerBound))  
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
