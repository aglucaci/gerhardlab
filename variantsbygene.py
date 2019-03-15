#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:43:25 2019

@author: alexander lucaci
"""

import sys
import plotly
import plotly.graph_objs as go
from plotly import tools
import numpy as np 
from scipy import stats

#Annotated VCF by VEP
fname = "PmhbHNI30Bkw48yH.vcf..txt"

gene_dict = {}


def plotly_things2(numbers, filename, desc1, desc2, desc3): #single historgram
    plotly.offline.init_notebook_mode(connected=False)
    plotly.offline.plot({"data": [go.Histogram(x=numbers, text=desc3)],"layout": go.Layout(title=desc2, xaxis={'title':desc1}, yaxis={'title':'Occurences'})},filename=filename+"_HISTOGRAM.html")

def plotly_things(chromos, numbers, filename, desc1, desc2): #single bar graph
    plotly.offline.init_notebook_mode(connected=False)
    plotly.offline.plot({"data": [go.Bar(x=chromos, y=numbers)],"layout": go.Layout(title=desc2, xaxis={'title':desc1}, yaxis={'title':'Occurences'})},filename=filename+"_HISTOGRAM.html")


def main(filename):
    global chr_dict
    with open(filename) as f:
        while True:
            line= f.readline().strip()
            if line == "": break
            if line[0] != "#": #skip header info
                #print(line.split()[7].split("|")[3])
                #print()
                

                try:
                    gene_dict[line.split()[7].split("|")[3]] += 1
                except:
                #    #print(line)                 
                    gene_dict[line.split()[7].split("|")[3]] = 1
    f.close()   
    
    
main(fname)

#print(gene_dict)

#take titles,
genes = [k for k in gene_dict]
numbers = [gene_dict[k] for k in gene_dict]

print(len(genes))
#print(chromos, numbers)

#plotly_things(genes, numbers, "WES_VCF_RUN_01",  "Gene", "Variants by Gene")
#plotly_things(TripleH_rates, "TH_RATES", "rate at which 3 nucleotides are changed instantly within a single codon", "Triple Mutation hit rate histogram", "Rate value")

# =============================================================================
# End of file
# =============================================================================
