#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: alexander lucaci
Gerhard LAB

provide gnomad AF filtering for VCF's from VEP

"""
# =============================================================================
#  Imports
# =============================================================================
import os, sys

# =============================================================================
#  Declares
# =============================================================================
fname = sys.argv[1]
with_gnomAD = 0
threshold = 0.00002
phenotype = []
gnomAD_AF_columninvcf = 0

# =============================================================================
#  Helper functions
# =============================================================================
def get_gnomadAF_loc(filename):
    global gnomAD_AF_columninvcf 
    with open(filename) as f:
        while True:
            line = f.readline().strip()
            
            if line == "": break

            if line[:6] == "##INFO": 
                if "VEP" in line:
                    INFO_Cols = line.split(",")[3].split("|")
                    for i, item in enumerate(INFO_Cols):
                    
                        if item == "gnomAD_AF": 
                            #print(i, [item])
                            gnomAD_AF_columninvcf = i

        f.close()

def apply_filter(data):
    global with_gnomAD, threshold, gnomAD_AF_columninvcf 
    msg = ""
    gnomAD_AF = data.split("\t")[7].split("|")[gnomAD_AF_columninvcf]
    if gnomAD_AF == "": 
        msg += data
    else:
        if float(gnomAD_AF) < threshold: 
            msg += data
            with_gnomAD += 1
    return msg
    
        
def read_vcf(filename):
    count, num_variants, header_count = 0, 0, 0
    
    with open(filename) as f:
        while True:
            line = f.readline().strip()
            if line == "": break

            if line[0] != "#":
                variant = apply_filter(line)

                if variant != "": 
                    #REGULAR VARIANTS 
                    print(variant)
                    num_variants += 1
            else:
                print(line)
                header_count += 1
            count += 1            
        f.close()
    return num_variants, count - header_count

# =============================================================================
#  MAIN          
# =============================================================================
#preprocessing to determing where to search for the gnomad_AF data
get_gnomadAF_loc(sys.argv[1])

total_vars, count = read_vcf(fname)

print()
print("() Total number of variants returned:", total_vars)
print("() Processed:", count)
print("() with gnomAD AF:", with_gnomAD)
print("() Threshold:", threshold)

for item in phenotype:
    print(item.split(" ")[0].replace("chr", "") + "[chr]", item.split(" ")[1] + ":" + item.split(" ")[1] +"[chrpos37]", item)


# =============================================================================
#  END OF FILE        
# =============================================================================
    
    
    
