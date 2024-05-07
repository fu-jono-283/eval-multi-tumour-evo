import sys
import os
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

IUPAC = {
        'R' : 'AG',
        'Y' : 'CT',
        'S' : 'CG',
        'W' : 'AT',
        'K' : 'GT',
        'M' : 'AC',
        'A' : 'A',
        'G' : 'G',
        'C' : 'C',
        'T' : 'T',
        'N' : 'ACGT'
}

consensus_iupac = {
        ('C', 'T'):'Y',
        ('A', 'G'): 'R',
        ('A', 'T'): 'W',
        ('C', 'G'): 'S',
        ('G', 'T'): 'K',
        ('A', 'C'): 'M',
        ('A', 'G', 'T'): 'N',
        ('A', 'C', 'G'): 'N',
        ('A', 'C', 'T'): 'N',
        ('C', 'G', 'T'): 'N',
        ('A', 'C', 'G', 'T'): 'N'
}

def consen(alig, ids, cutoff):
    s_bases = ('A', 'C', 'G', 'T')
    d_bases = ('R', 'Y', 'S', 'W', 'K', 'M', 'N')
    bases = s_bases + d_bases
    con_seq = ""
    for i in range(alig.get_alignment_length()):
        col = alig[:, i]
        col = col.upper()
        base_count = {b: col.count(b) for b in bases}
        genotypes = sorted(base_count.items(), key=lambda x: -x[1])
        # genotypes = [b for b in genotypes if b[1] > 0]
        genotypes = [b for b in genotypes if b[1] >= len(col) * cutoff]
        # genotypes = [e for e in genotypes if not (genotypes[0][1] !=  e[1] and e[1] <= 1)]
        if len(genotypes) == 0:
            con_seq += 'N'
        else:
            if len(genotypes) <= 1:
                con_seq += genotypes[0][0]
            else:
                amb_base = [IUPAC[c[0]] for c in genotypes]
                amb_base_list = []
                for d in amb_base:
                    amb_base_list.extend(list(d))
                amb_base_set = sorted(set(amb_base_list))
                con_seq += consensus_iupac[tuple(amb_base_set)]

    con_rec = SeqRecord(Seq(con_seq), id=ids, description='')
    return con_rec


if __name__ == "__main__":
    cutoff = float(sys.argv[1])  # Get the cutoff value from the command-line argument
    input_file = sys.argv[2]

    ids = []
    ids_con = []

    for seq_re in SeqIO.parse(input_file, "phylip"):
        ids.append(seq_re.id.split('_')[0])
    ids = sorted(set(ids))

    for key in ids:
        alignments = [seq for seq in SeqIO.parse(input_file, "phylip") if seq.id.split('_')[0] == key]
        alignments = MultipleSeqAlignment(alignments)
        ids_con.append(consen(alignments, key, cutoff))

    output_file = f"{input_file.split('.')[0]}_consensus_{cutoff:.2f}.fasta"  # Customize the output filename
    SeqIO.write(ids_con, output_file, "fasta-2line")