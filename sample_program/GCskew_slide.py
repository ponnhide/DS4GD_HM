from Bio import SeqIO 
from Bio import SeqUtils 
if __name__ == "__main__": 
    size = 10000
    for fna in SeqIO.parse("../sample_genome/sample2.fasta","fasta"): 
        results = []
        seq = fna.seq
        results = SeqUtils.GC_skew(seq,size) 
        with open("{}_GCskew.txt".format(fna.id),"w") as o:
            o.write(fna.id + "," + str(size) + "\n")
            for i in range(len(results)):
                o.write(str(i*size+1) + "," + str(results[i]) + "\n")
