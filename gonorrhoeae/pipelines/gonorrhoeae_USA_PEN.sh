#!/usr/bin/env bash

######################################################################################
##
## FIRST SECTION:
## Here we define the requirements we need to run our job: CPUs, RAM, run time, ...
##
######################################################################################

#BSUB -J gon_USA_pen # job name
#BSUB -n 10          # number of tasksi/CPUs in job
#BSUB -q standard   # queue
#BSUB -e error.%J   # error file name in which %J is replaced by the job ID
#BSUB -o output.%J  # output file name in which %J is replaced by the job ID

######################################################################################
##
## SECOND SECTION:
## Here we define the modules we need to be loaded in order to run our job.
## In this case, we need to load the R module, ... but first we make sure
## we start off with a sane, blank environment so we purge all modules
##
######################################################################################

module purge
module load r/4.2.2

######################################################################################
##
## THIRD SECTION:
## Here we define the executable/binary/script we want LSF to run, our real data
## crunching section. Remember, in our case: we want to run an R script that
## calculates factorial of 8 and our R script is called "r.script".
##
## As you can see, we specify the FULL path to the script and we redirect the output
## to a file.
##
######################################################################################

R < /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines/gonorrhoeae_USA_PEN.R --save >> /nfs/research/jlees/jacqueline/gwas_datasets/gonorrhoeae/pipelines/gon_USA_pen_results
