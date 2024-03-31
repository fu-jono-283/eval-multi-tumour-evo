import vcf
import sys
import os
from collections import namedtuple

vcf_reader = vcf.Reader(open(sys.argv[1], 'r'))

vcf_writer = vcf.Writer(open(sys.argv[1]+"_consensus.vcf", 'w'), vcf_reader)

PL_GT=['0/0','0/1','1/1','0/2','1/2','2/2','0/3','1/3','2/3','3/3']

def sum_lists(*args):
    return list(map(sum, zip(*args)))

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
        max_PL_value = max(sum_lists(*PL_values))
        combined_sample = vcf.model._Call(sample=key,site=combined_record,data=combined_record)
        combined_CallData = namedtuple('CallData', ["GT","DP","RC","G10","PL"])
        combined_sample.data = combined_CallData(GT=PL_GT[sum_lists(*PL_values).index(max_PL_value)], DP=DP_values, RC=sum_lists(*RC_values), G10=sum_lists(*G10_values), PL=sum_lists(*PL_values))
        con.append(combined_sample)
    combined_record.samples=con
    combined_record.FORMAT="GT:DP:RC:G10:PL"
    #remove record if  multi allele bigger than 2
    if len(sys.argv) == 2:
        vcf_writer.write_record(combined_record)
    elif sys.argv[2] == "snp2" and len(set("".join(genotypes)))==2:
        vcf_writer.write_record(combined_record)

vcf_writer.close()

with open(sys.argv[1]+"_consensus.vcf", 'r') as vcf_file:
    update_header=[]
    for j in vcf_file.readlines():
        if j.startswith('#CHROM'):
            header=j.split()
            index=header.index("FORMAT")
            header[index + 1:]=ids
            update_header.append('\t'.join(map(str,header)) + '\n')
        else:
            update_header.append(j)
with open(sys.argv[1]+"_consensus.vcf", 'w') as vcf_file:
    vcf_file.writelines(update_header)

