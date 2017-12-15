# TopicModels
# Abstract
This tools can estimate parameters(including hyper parameters) of some topic models.  
# Installation
    ./configure
    make
# Requirements
You should install boost library for using this tools.  

If you also use scripts, you should install miniconda.  
After installation of it, you should execute following commands.  
`conda install numpy scipy pandas`  
# Usage
## Latent Dirichlet Allocation
 LDA [BOW file] [-options]  
  
Options:  
  - -h [ --help ]:             show help  
  - -k [ --nmtp ] arg (=10):   number of topics  
  - -f [ --nmsh ] arg (=5):    number of factors with high probability to show  
  - -o [ --otpt ] arg (=./):   directory name for output  
  
# Input
  Input file is required to be frequency matrix format. See data/test.csv.  
# Output
  This program outputs following:  
## LDA
  1)topic distribution(theta.csv)  
  2)word distribution(phi.csv)  
  3)hyper parameters(alpha.csv, beta.csv)  
  4)word list(wordList)  
  ( 5)vatiational lower bound(variationalLowerBound))  
# Scripts
  You can also use python scripts. These scripts calculate some values of each number of topics.  
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
  
deriveldavlb.py [-h] [-s KSTEP] BOW_filename k_min k_max  
  
### Maximized variational lower bound
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
These scripts output results to /output/(time stamp) directory. (time stamp) means date information of MMDD_DAY_hh_mm_ss format.  
# License
This software is released under the MIT License, see LICENSE.txt.  
