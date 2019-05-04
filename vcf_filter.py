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
consequence = ["synonymous_variant", "intron_variant", "non_coding_transcript_exon_variant", "upstream_gene_variant", "downstream_gene_variant"]
fname = "5enxcsqdMxN4a9CH.vcf"

# =============================================================================
# Helper functions
# =============================================================================
def apply_filter(data):
    for filter_item in consequence:
        if filter_item in data:
            return ""
        
        #OR
        
        try:
            x = float(data.split("\t")[7].split("|")[42])
            if x == "": print("NADA")
            if x > 0.01:
                #print([x])
                return ""
        except:
            pass
        
    return data



def read_vcf(filename):
    count = 0
    header = ""
    variants = ""
    numvariants = 0
    
    with open(filename) as f:
        while True:
            line = f.readline().strip()
            if line == "": break
            if count == 10000000000: break #start at 399
            
            header += line + "\n"
            
            if line[0] != "#":
                x = apply_filter(line)
                #if x == "": continue
                print(x)
                
                numvariants += 1
                #print(apply_filter(line).split("\t")[7].split("|")[42])
                #print(line)
                
            count += 1
            
        f.close()
    return numvariants, count
# =============================================================================
#  MAIN          
# =============================================================================
            
total_vars, count = read_vcf(fname)

print("Total number of variants returned:", total_vars, count)