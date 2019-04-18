"""
 Creates all tables and figures for Chapter 1 of my dissertation
"""

# Requrements:
#     numdifftools
#     numpy
#     matplotlib
#     scipy

from __main__ import *
import subprocess

#Rscript_path = "C:/Program Files/R/R-3.4.3/bin/x64/Rscript"

print('First create tables and graphs in Python')
import CreateTable6

print('Next run the R script that creates figures 1 and 2')
r_status = subprocess.call([Rscript_path, "Rcode/time_agg_random_walk_graph.R"], shell=True)
if r_status!=0:
    print('R code could not run. Check the the path for Rscript')
    
    
