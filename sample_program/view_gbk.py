#get_CDS.py
from Bio import SeqIO 
def get_CDS_info(gbk):
    for feat in gbk.features:
        start = feat.location.parts[0].start.position
        end   = feat.location.parts[-1].end.position
        print(start+1,end,feat.strand,feat.type)
        for key in feat.qualifiers:
            print(key,feat.qualifiers[key][0])
        print() 

if __name__ == "__main__":
    for gbk in SeqIO.parse("../sample_genome/sample.gbk","genbank"):
    	get_CDS_info(gbk)
