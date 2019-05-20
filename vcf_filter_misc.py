#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  4 15:26:48 2019

@author: alexander lucaci
Gerhard LAB

provide gnomad AF filtering for VCF's from VEP

"""

"""

##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|TSL|APPRIS|ENSP|SWISSPROT|TREMBL|UNIPARC|REFSEQ_MATCH|GIVEN_REF|USED_REF|BAM_EDIT|SIFT|PolyPhen|HGVS_OFFSET|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|gnomAD_OTH_AF|gnomAD_SAS_AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">


"""
import os, sys

# =============================================================================
# Declares
# =============================================================================
#consequence = ["synonymous_variant", "intron_variant", "non_coding_transcript_exon_variant", "upstream_gene_variant", "downstream_gene_variant"]
#fname = "jVTFhw6xfxgKfarh.Consequence_ne_downstream_gene_variant_and_Consequence_ne_upstream_gene_variant_and_Consequence_ne_intron_variant_and_Consequence_ne_synonymous_variant_and_Consequence_ne_non_coding_transcript_exon_variant_and_Consequence_n.vcf"
#fname = "MoreStringent_jVTFhw6xfxgKfarh.Consequence_ne_downstream_gene_variant_and_Consequence_ne_upstream_gene_variant_and_Consequence_ne_intron_variant_and_Consequence_ne_synonymous_variant_and_Consequence_ne_non_coding_transcript_exon_variant_and_Consequen.vcf"
fname = "jVTFhw6xfxgKfarh.Consequence_is_missense_variant.vcf"
#fname = "jVTFhw6xfxgKfarh.Consequence_is_missense_variant_and_clinvar_clnsig_is_5.vcf"
with_gnomAD = 0
threshold = 0.00002

phenotype = []

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
            #if count == 500: break  #starts at 333
            
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
                    #REGULAR VARIANTS 
                    """print(x) #FOR VCF MODE """
                    
                    
                    
                    #print(x.split("\t")[0].replace("chr", "") + "[chr] AND", x.split("\t")[1] + ":" + x.split("\t")[1])
                    #https://www.ncbi.nlm.nih.gov/clinvar/?term=6%5Bchr%5D+AND+9708032%3A9708032
                    #print(x.split("\t"))
                    
                    #if x.split("\t")[1] == str(152282024): 
                    #    print(x.split("\t")[7].split("|")[231]) # is clinvar
                    
                    
                    
                    #CLINVAR
                    """
                    if x.split("\t")[7].split("|")[231] != "":
                        #print([x.split("\t")[7].split("|")[231]])
                        print(x)
                    
                        #print()
                    
                        num_variants += 1
                        
                        #if "5" in x.split("\t")[7].split("|")[231]:
                        phenotype.append(x.split("\t")[0] + " " + x.split("\t")[1] + " " + x.split("\t")[3] + ">" + x.split("\t")[3] + "|" + x.split("\t")[7].split("|")[231])
                    """ 
                        
                        
                        
                    num_variants += 1
                    
            count += 1
            
        f.close()
    return num_variants, count - 333
# =============================================================================
#  MAIN          
# =============================================================================
            
total_vars, count = read_vcf(fname)

print("() Total number of variants returned:", total_vars)
print("() Processed:", count)
print("() with gnomAD AF:", with_gnomAD)
print("() Threshold:", threshold)

for item in phenotype:
    print(item.split(" ")[0].replace("chr", "") + "[chr]", item.split(" ")[1] + ":" + item.split(" ")[1] +"[chrpos37]", item)
    
    
    
