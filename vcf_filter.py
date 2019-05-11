#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  4 15:26:48 2019

@author: iamtokr
"""


import os, sys

# =============================================================================
# Declares
# =============================================================================
#consequence = ["synonymous_variant", "intron_variant", "non_coding_transcript_exon_variant", "upstream_gene_variant", "downstream_gene_variant"]
#fname = "jVTFhw6xfxgKfarh.Consequence_ne_downstream_gene_variant_and_Consequence_ne_upstream_gene_variant_and_Consequence_ne_intron_variant_and_Consequence_ne_synonymous_variant_and_Consequence_ne_non_coding_transcript_exon_variant_and_Consequence_n.vcf"
#fname = "MoreStringent_jVTFhw6xfxgKfarh.Consequence_ne_downstream_gene_variant_and_Consequence_ne_upstream_gene_variant_and_Consequence_ne_intron_variant_and_Consequence_ne_synonymous_variant_and_Consequence_ne_non_coding_transcript_exon_variant_and_Consequen.vcf"
fname = "jVTFhw6xfxgKfarh.Consequence_is_missense_variant.vcf"
with_gnomAD = 0
threshold = 0.00002

# =============================================================================
# Helper functions
# =============================================================================
def apply_filter(data):
    global with_gnomAD, threshold
    msg = ""
    #return data.split("\t")[7].split("|")[43]
    gnomAD_AF = data.split("\t")[7].split("|")[43]
    #print(type(gnomAD_AF))
    #try:
        #print(gnomad_AF.isdigit())
    #    return float(gnomAD_AF)
    #except:
    #    return(data)
    
    if gnomAD_AF == "": 
        #return data
        msg += data
    else:
        #return(float(gnomAD_AF))
        #if float(gnomAD_AF) < .5: return data
        if float(gnomAD_AF) < threshold: 
            msg += data
            with_gnomAD += 1
    return msg
    
        
def read_vcf(filename):
    count = 0
    #header = ""
    num_variants = 0
    
    with open(filename) as f:
        while True:
            line = f.readline().strip()
            if line == "": break
            #if count == 333: break  #starts at 333
            
            if line[0] == "#": 
                #header += line + "\n"
                print(line)
                
            if line[0] != "#":

                #print(line)
                x = apply_filter(line)
                #try:
                #    print([x], len(x))
                #except:
                #    print([x])
                
                if x != "": 
                    print(x)
                
                    #print()
                    num_variants += 1

            count += 1
            
        f.close()
    return num_variants, count - 333
# =============================================================================
#  MAIN          
# =============================================================================
            
total_vars, count = read_vcf(fname)

#print("() Total number of variants returned:", total_vars)
#print("() Processed:", count)
#print("() with gnomAD AF:", with_gnomAD)
#print("() Threshold:", threshold)