from Bio import SeqIO 
from Bio import SeqUtils 
if __name__ == "__main__": 
    for fna in SeqIO.parse("../sample_genome/sample1.fasta","fasta"): 
        print(fna.id)
        print(fna.seq) 
        print(fna.seq.reverse_complement()) 
        print(SeqUtils.GC(fna.seq)) 
        print(SeqUtils.GC_skew(fna.seq)[0])
