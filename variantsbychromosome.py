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


fname = "SQ8992_L001_R1_R2_hg38_sorted.mpileup_variants.vcf"

chr_dict = {}

for n in range(1,23):
    #print("chr"+str(n))
    chr_dict["chr"+str(n)] = 0
chr_dict["chrX"] = 0
chr_dict["chrY"] = 0

#print(chr_dict)

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
                #print(line.split()[0])
                try:
                    chr_dict[line.split()[0]] += 1
                except:
                    #print(line)                 
                    chr_dict[line.split()[0]] = 1  
    f.close()   
main(fname)

#print(chr_dict)

#take titles,
chromos = [k for k in chr_dict]
numbers = [chr_dict[k] for k in chr_dict]

#print(chromos, numbers)

plotly_things(chromos, numbers, "WES_VCF_RUN_01",  "Chromosome", "Variants by chromosome")



#plotly_things(TripleH_rates, "TH_RATES", "rate at which 3 nucleotides are changed instantly within a single codon", "Triple Mutation hit rate histogram", "Rate value")


