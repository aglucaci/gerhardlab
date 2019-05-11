#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 15 11:43:25 2019

@author: alexander lucaci
"""
# =============================================================================
# Imports
# =============================================================================
import sys
import plotly
import plotly.graph_objs as go
from plotly import tools
import numpy as np 
from scipy import stats
from prettytable import PrettyTable
#http://zetcode.com/python/prettytable/

# =============================================================================
# Declares
# =============================================================================
#Annotated VCF by VEP
#path = "/Users/phylo/Documents/Gerhard Lab/RQ534361/"
#fname = path + "zFSM84RbEYaQM8ck.gnomAD_AF_lt_0.01_and_Consequence_ne_downstream_gene_variant_and_Consequence_ne_upstream_gene_variant_and_Consequence_ne_intron_variant_and_Consequence_ne_synonymous_variant_and_Consequence_ne_non_coding_transcript_exon_var.vcf"
#path = "/Users/phylo/Documents/Gerhard\ Lab/RQ534361/VEP Analysis 05042019/"

#fname = "jVTFhw6xfxgKfarh.Consequence_ne_downstream_gene_variant_and_Consequence_ne_upstream_gene_variant_and_Consequence_ne_intron_variant_and_Consequence_ne_synonymous_variant_and_Consequence_ne_non_coding_transcript_exon_variant_and_Consequence_n.vcf"
#fname = "gnomAD_AF_Filtered_jVTFhw6xfxgKfarh_Consequences.vcf"
fname = "Removed_chrY_gnomADAF_flt_Alex.vcf"
gene_dict = {}

# =============================================================================
# Helper functions
# =============================================================================
def plotly_things2(numbers, filename, desc1, desc2, desc3): #single historgram
    plotly.offline.init_notebook_mode(connected=False)
    plotly.offline.plot({"data": [go.Histogram(x=numbers, text=desc3)],"layout": go.Layout(title=desc2, xaxis={'title':desc1}, yaxis={'title':'Occurences'})},filename=filename+"_HISTOGRAM.html")

def plotly_things(chromos, numbers, filename, desc1, desc2): #single bar graph
    plotly.offline.init_notebook_mode(connected=False)
    plotly.offline.plot({"data": [go.Bar(x=chromos, y=numbers)],"layout": go.Layout(title=desc2, xaxis={'title':desc1}, yaxis={'title':'Occurences'})},filename=filename+"_HISTOGRAM.html")

def plotly_things_vert(chromos, numbers, filename, desc1, desc2): #single bar graph
    plotly.offline.init_notebook_mode(connected=False)
    plotly.offline.plot({"data": [go.Bar(y=chromos, x=numbers)],"layout": go.Layout(title=desc2, xaxis={'title':"Occurences"}, yaxis={'title':desc1})},filename=filename+"_HISTOGRAM.html")

# =============================================================================
# main sub
# =============================================================================
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
    
# =============================================================================
# main prog. starts here
# =============================================================================
main(fname)

#print(gene_dict)
#take titles,
#genes = [k for k in gene_dict]
#numbers = [gene_dict[k] for k in gene_dict]
#print(len(genes))
#plotly_things(genes, numbers, "WES_VCF_RUN_01",  "Gene", "Variants by Gene")
#plotly_things(TripleH_rates, "TH_RATES", "rate at which 3 nucleotides are changed instantly within a single codon", "Triple Mutation hit rate histogram", "Rate value")

# --- Sort dict for histo ----- #
#print(sorted(gene_dict, key=gene_dict.get))

d = gene_dict
genes, numbers = [], []
#is this accurate?
for w in sorted(d, key=d.get, reverse=True):
  #print(w, d[w])
  genes.append(w)
  numbers.append(d[w])
  
#genes[0] = "UNKNOWN"

#genes.pop(0)
#numbers.pop(0)

print("Filename:", fname)
print("Number of genes:", len(genes), ",Number of variants:", sum(numbers))
plotly_things(genes, numbers, "WES_VCF_RUN_01",  "Gene", "Variants by Gene, (Gene n=" + str(len(genes)) + ") (Variants n=" + str(sum(numbers)) + ")")
#print(*genes[:25], sep="\n")


x = PrettyTable()
x.field_names = ["#", "Gene", "Number of variants"]
                 
for i, gene in enumerate(genes[:25]):
    print(i+1, "GENE:", gene, ", Variants:", numbers[i])
    x.add_row([i+1, gene, numbers[i]])
print()
print(x)
    
x = PrettyTable()
x.field_names = ["#", "Gene", "Number of variants"]
  
                
for gene in genes: print(gene)
print(",".join(genes))
print(" OR ".join(genes))


# =============================================================================
# End of file
# =============================================================================
