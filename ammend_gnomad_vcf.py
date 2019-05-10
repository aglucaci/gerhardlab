# -*- coding: utf-8 -*-

"""

@author: alexander lucaci

Corrects the naming convention for gnomad VCF

1 -> chr1 
so that it agrees with hg19 reference.
"""
import os, sys

# =============================================================================
# Declares
# =============================================================================
fname = "gnomad.exomes.r2.1.1.sites.vcf"

# =============================================================================
# Helper functions
# =============================================================================



def read_vcf(filename):
    count = 0
    header = ""
    variants = ""
    numvariants = 0
    
    
    with open(filename) as f:
        while True:
            line = f.readline().strip()
            
            if line == "": break
            
            #if count == 2: break
        
            
        
            if line[0] != "#":
                holder = []
                #print("chr"+ line.split("\t")[0])
                holder.append("chr"+ line.split("\t")[0])
                holder += line.split("\t")[1:]
                
                print("\t".join(holder))
                count += 1
            else:
                print(line)
            
        f.close()

# =============================================================================
#  MAIN          
# =============================================================================
            
read_vcf(fname)

