
from Bio import SeqIO





inpath = '/home/luotao/c_elegans/WGS/ref/Glycine_max.Glycine_max_v2.1.dna.toplevel.fa'

outpath = '/home/luotao/c_elegans/WGS/ref/Glycine_max.Glycine_max_v2.1.dna.toplevel.intervals'




record = SeqIO.parse(inpath, "fasta")


with open(outpath,'w') as f0:
    for seq_record in record:
        f0.write(f'{seq_record.description}:1-{str(len(seq_record))}\n')