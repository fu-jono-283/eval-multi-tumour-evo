    # This script is used to produce a consensus sequence fasta file from the CellCoal haplotype file. The cutoff value can be adjusted. 

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

    def consen(alig, ids):
        s_bases = ('A', 'C', 'G', 'T')
        d_bases = ('R', 'Y', 'S', 'W', 'K', 'M', 'N')
        bases = s_bases + d_bases
        con_seq = ""
        cutoff = 0.35
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


    ids = []
    ids_con = []

    for seq_re in SeqIO.parse(sys.argv[1], "phylip"):
        ids.append(seq_re.id.split('_')[0])
        ids=sorted(set(ids))

    for key in ids:
        alignments=[seq for seq in SeqIO.parse(sys.argv[1], "phylip") if seq.id.split('_')[0]==key]
        alignments=MultipleSeqAlignment(alignments)
        ids_con.append(consen(alignments,key))

    SeqIO.write(ids_con, sys.argv[1]+"_consensus.fasta", "fasta-2line")