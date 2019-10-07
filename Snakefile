__author__ = "Aaron Doe"
__email__ = "doe@ohsu.edu"
__license__ = "MIT"

"""Computation Hub omic data processing pipeline for SE Bulk RNA-sequencing"""

import datetime
import sys
import os
import pandas as pd
import json

timestamp = ('{:%Y-%m-%d_%H:%M:%S}'.format(datetime.datetime.now()))

with open('cluster.json') as json_file:
     json_dict = json.load(json_file)

rule_dirs = list(json_dict.keys())

for rule in rule_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'logs',rule)):
       log_out = os.path.join(os.getcwd(), 'logs', rule)
       os.makedirs(log_out)
       print(log_out)


result_dirs = ['diffexp','tables']
for rule in result_dirs:
    if not os.path.exists(os.path.join(os.getcwd(),'results',rule)):
       log_out = os.path.join(os.getcwd(), 'results', rule)
       os.makedirs(log_out)
       print(log_out)


def message(mes):
    sys.stderr.write("|--- " + mes + "\n")


rule all:
    input:
       "plots/cor_acc/acc_rna_correlations.pdf",
       "plots/cor_met/met_rna_correlations.pdf",
       "plots/cor_accmetrna/accmetrna_correlations.scatter.pdf",
       "plots/cor_accmetrna/accmetrna_correlations.tsv",
       "plots/cor_accmetrna/accmetrna_correlations.boxplot.pdf"
        

include: "rules/correlations.smk"
