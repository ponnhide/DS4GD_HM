def read_fasta(fname):
    seq_dict = {} 
    with open(fname,"r") as f:
        for line in f:
            if line[0] == ">":
                seq_id = line.rstrip()[1:]
                seq_dict[seq_id] = ""
            else:
                seq_dict[seq_id] += line.rstrip().upper()#念のため大文字化 
    return seq_dict

def calc_GC(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) 
    
def calc_GC_skew(seq):
    G_num = seq.count("G")
    C_num = seq.count("C")
    return (G_num - C_num)/ (G_num + C_num)
    
if __name__ == "__main__": 
    seq_dict = read_fasta("../sample_genome/sample1.fasta")
    for seq_id in seq_dict.keys(): 
        seq = seq_dict[seq_id]
        print(seq)
        print(seq.translate(str.maketrans("ATGC","TACG"))[::-1]) 
        print(calc_GC(seq))
        print(calc_GC_skew(seq))
