from Bio import SeqIO 
from Bio import SeqUtils 
if __name__ == "__main__": 
    size = 10000
    for fna in SeqIO.parse("../sample_genome/sample2.fasta","fasta"): 
        results = []
        seq = fna.seq
        for i in range(0,len(seq),size):
            GC_ratio = SeqUtils.GC(seq[i:i+size])
            results.append(GC_ratio)
        
        with open("{}_GCratio.txt".format(fna.id),"w") as o:
            o.write(fna.id + "," + str(size) + "\n")
            for i in range(len(results)):
                o.write(str(i*size+1) + "," + str(results[i]) + "\n")
