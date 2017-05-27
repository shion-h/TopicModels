# TopicModels
# Abstract
This tools can estimate parameters(includes hyper parameters) of some topic models.  
# Environment
# Installation
You should install boost library for using this tools.  
And then, execute make command at root directory of tools. After that, execution files exist in /bin/ directory.  

If you also use scripts, you should install miniconda.  
After installation of it, you should execute following commands.  
    conda install numpy, scipy, pandas  

# Usage
## Latent Dirichlet Allocation
 LDA [BOW file] [-options]  
  
Options:  
  -h [ --help ]             show help  
  -k [ --nmtp ] arg (=10)   number of topics  
  -s [ --nmit ] arg (=1000) number of iterations  
  -b [ --bnin ] arg         burn in period  
  -i [ --intv ] arg (=10)   sampling interval  
  -l [ --lrna ] arg (=1)    learning algorythm(0:Gibbs sampling 1:Collapsed gibbs sampling 2:Variational Bayes)  
  -f [ --nmsh ] arg (=5)    number of factors with high probability to show  
  -o [ --otpt ] arg (=./)   directory name for output  
  
## Author-Topic Model
 ATM [BOW file] [Author file] [-options]  
  
Options:  
  -h [ --help ]             show help  
  -k [ --nmtp ] arg (=10)   number of topics  
  -s [ --nmit ] arg (=1000) number of iterations  
  -b [ --bnin ] arg         burn in period  
  -i [ --intv ] arg (=10)   sampling interval  
  -f [ --nmsh ] arg (=5)    number of factors with high probability to show  
  -o [ --otpt ] arg (=./)   directory name for output  
# Input
  Input file is required to be frequency matrix format. See data/test.csv.  
# Output
  This program outputs following:  
## LDA
  1)topic distribution(theta.csv) 2)word distribution(phi.csv) 3)hyper parameters(alpha.csv, beta.csv) 4)word list(wordList) ( 5)vatiational lower bound(variationalLowerBound))  
## ATM
  1)topic distribution(theta.csv) 2)word distribution(phi.csv) 3)hyper parameters(alpha.csv, beta.csv) 4)word list(wordList) 5)author list (authorList)  
# Scripts
  You can also use python scripts. These scripts calculate some values of each number of topics.  
## Usage
calcldappl.py [-h] [-s KSTEP] [-r RATE] [-i ITERATION]  
                     BOW_filename k_min k_max  
  
positional arguments:  
  BOW_filename          Bag Of Words filename  
  k_min                 minimum number of topics for calculation of perplexity  
  k_max                 max number of topics for calculation of perplexity  
  
optional arguments:  
  -h, --help            show this help message and exit  
  -s KSTEP, --kstep KSTEP step number of topics for calculation of perplexity  
  -r RATE, --rate RATE  rate of train data  
  -i ITERATION, --iteration ITERATION number of iterations(>10)  
  
calcldawbic.py [-h] [-s KSTEP] [-i ITERATION] BOW_filename k_min k_max  
  
positional arguments:  
  BOW_filename          Bag Of Words filename  
  k_min                 minimum number of topics for calculation of wbic  
  k_max                 max number of topics for calculation of wbic  
  
optional arguments:  
  -h, --help            show this help message and exit  
  -s KSTEP, --kstep KSTEP step number of topics for calculation of wbic  
  -i ITERATION, --iteration ITERATION number of iterations(>10)  
  
calcldavlb.py [-h] [-s KSTEP] BOW_filename k_min k_max  
  
positional arguments:  
  BOW_filename          Bag Of Words filename  
  k_min                 minimum number of topics for calculation of maximized variational lower bound  
  k_max                 max number of topics for calculation of maximized variational lower bound  
  
optional arguments:  
  -h, --help            show this help message and exit  
  -s KSTEP, --kstep KSTEP step number of topics for calculation of maximized variational lower bound  
## Output
These scripts output results to /output/(time stamp) directory. (time stamp) means date information of MMDD_DAY_hh_mm_ss format.  
# License
This software is released under the MIT License, see LICENSE.txt.  

