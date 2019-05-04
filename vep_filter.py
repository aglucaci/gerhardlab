#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  4 16:38:17 2019

@author: iamtokr
"""

#fname = "5enxcsqdMxN4a9CH.Consequence_ne_intron_variant_and_Consequence_ne_downstream_gene_variant_and_Consequence_ne_synonymous_variant_and_Consequence_ne_non_coding_transcript_exon_variant_and_Consequence_ne_upstream_gene_variant.vep.txt"

fname = "hkiQ0C5ETbqbXztT.Consequence_ne_downstream_gene_variant_and_Consequence_ne_upstream_gene_variant_and_Consequence_ne_intron_variant_and_Consequence_ne_synonymous_variant_and_Consequence_ne_non_coding_transcript_exon_variant_and_Consequence_n.txt"
gene_dict = {}

def filtered_genedict(data):
    #Genedict
    if "SYMBOL=" in item.upper():
        gene = item.upper().replace("SYMBOL=", "")
        try:
            gene_dict[gene] += 1
        except:
            gene_dict[gene] = 1
                    
        if "gnomAD_AF=" in item: 
            x = float(item.replace("gnomAD_AF=", ""))

def read_vep(filename):
    count = 0
    header = ""
    variants = ""
    num_variants = 0
    
    with open(filename) as f:
        while True:
            line = f.readline().strip()
            #count += 1
            
            
            if line == "": break
            if count == 2000000: break #start at 399
            
            if line[0] == "#": 
                header += line + "\n"
                print(line)
                
            if line[0] != "#":
                #print(line.split("\t")[13], "\n")
                for item  in line.split("\t")[13].split(";"):
                        if "gnomAD_AF=" in item: 
                            x = float(item.replace("gnomAD_AF=", ""))
                        
                            if x < 0.01: 
                                #print(line)
                                num_variants += 1
                                #break
                            
                                for item  in line.split("\t")[13].split(";"):
                                    if "SYMBOL=" in item.upper():
                                        gene = item.upper().replace("SYMBOL=", "")
                                        try:
                                            gene_dict[gene] += 1
                                        except:
                                            gene_dict[gene] = 1
                            
                            
                            
                            
                            
                            
                        
            
            
            
        f.close()
        
    return num_variants, count

# =============================================================================
#  MAIN          
# =============================================================================

num_variants, count = read_vep(fname)

print(num_variants, count)


d = gene_dict
genes, numbers = [], []
#is this accurate?
for w in sorted(d, key=d.get, reverse=True):
  #print(w, d[w])
  genes.append(w)
  numbers.append(d[w])
  
genes[0] = "UNKNOWN"
#genes.pop(0)
#numbers.pop(0)

print("Number of genes:", len(genes), ",Number of variants:", sum(numbers))
#plotly_things(genes, numbers, "WES_VCF_RUN_01",  "Gene", "Variants by Gene, (Gene n=" + str(len(genes)) + ") (Variants n=" + str(sum(numbers)) + ")")
#print(*genes[:25], sep="\n")

msg = ""

for i, gene in enumerate(genes[:25]):
    print(i+1, gene)
    msg += gene + " OR "
print(msg)

for i, gene in enumerate(genes[:25]):
    print(gene)