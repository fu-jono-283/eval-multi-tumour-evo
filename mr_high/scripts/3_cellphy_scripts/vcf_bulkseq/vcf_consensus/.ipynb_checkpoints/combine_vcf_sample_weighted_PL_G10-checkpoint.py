#####This script is for combining multi samples in a VCF file, and use the normalized value of PL to detect the consensus genotype.
#####Normolized process: for each snp site, for each combined species, sum the PL value for each sample, then use (1/this value) times each PL value of this sample; then sum all the combined samples together by using this formula,then divide by the total value of all the combined samples on this site.   
import vcf
import sys
import os
import math
from collections import namedtuple

vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))

if sys.argv[1].endswith(".vcf"):
    output_file = sys.argv[1].replace(".vcf","_consensus.vcf.NEW_ER")
else:
    output_file = sys.argv[1] + "_consensus.vcf.NEW_ER" 

vcf_writer = vcf.Writer(open(output_file, 'w'), vcf_reader)

PL_GT=['0/0','0/1','1/1','0/2','1/2','2/2','0/3','1/3','2/3','3/3']

def sum_lists(*args):
    return list(map(sum, zip(*args)))

def weighted_PL(PL_list):
    #has_sublists = any(isinstance(item, list) for item in PL_list)
    if len(PL_list) ==1:
        return PL_list[0]
    sublist_sums = [sum(xx for xx in sublist if xx!= -math.inf and xx is not None and xx!=-2147483648) for sublist in PL_list]
    total_sum = sum(1 / total for total in sublist_sums if total != 0)
    indiv_res = [[(1 / total) * value if total != 0 else 0 for value in sublist] for sublist, total in zip(PL_list, sublist_sums)]
    #indiv_res = [[(1/sum(sublist)) * value if sum(sublist) !=0 and value != math.inf and value != "." else 0 for value in sublist] for sublist in PL_list]
    norm_res = [round(x/total_sum,2) for x in sum_lists(*indiv_res)]
    return norm_res

ids = []
for i in vcf_reader.samples:
    ids.append(i.split('_')[0])
    ids=sorted(set(ids))

for record in vcf_reader:
    combined_record = record
    con=[]
    genotypes = []
    for key in ids:
        DP_values=0
        PL_values=[]
        RC_values=[]
        G10_values=[]
        for sample in record.samples:
            sep = ('/' if not sample.phased else '|')
            if sample.gt_bases != None:
                genotypes.append(''.join(sample.gt_bases.replace(sep,'')))
            if sample.sample.split('_')[0]==key:
                DP_values+=sample.data.DP
                PL_values.append(sample.data.PL)
                RC_values.append(sample.data.RC)
                G10_values.append(sample.data.G10)
        max_PL_value = max(weighted_PL(PL_values))
        combined_sample = vcf.model._Call(sample=key,site=combined_record,data=combined_record)
        combined_CallData = namedtuple('CallData', ["GT","DP","RC","G10","PL"])
        combined_sample.data = combined_CallData(GT=PL_GT[weighted_PL(PL_values).index(max_PL_value)], DP=DP_values, RC=sum_lists(*RC_values), G10=weighted_PL(G10_values), PL=[round(n) for n in weighted_PL(PL_values)])
        con.append(combined_sample)
    combined_record.samples=con
    combined_record.FORMAT="GT:DP:RC:G10:PL"
    ##remove record if  multi allel bigger than 2
    if len(sys.argv) == 2:
        vcf_writer.write_record(combined_record)
    elif sys.argv[2] == "snp2" and len(set("".join(genotypes)))==2:
        vcf_writer.write_record(combined_record)

vcf_writer.close()

with open(output_file, 'r') as vcf_file:
    update_header=[]
    for j in vcf_file.readlines():
        if j.startswith('#CHROM'):
            header=j.split()
            index=header.index("FORMAT")
            header[index + 1:]=ids
            update_header.append('\t'.join(map(str,header)) + '\n')
        else:
            update_header.append(j)
with open(output_file, 'w') as vcf_file:
    vcf_file.writelines(update_header)

