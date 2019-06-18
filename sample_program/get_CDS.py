#get_CDS2.py
from Bio import SeqIO 
def get_CDS_info(gbk):
    results = []
    for feat in gbk.features:
        start  = feat.location.parts[0].start.position #Featureの開始位置
        end    = feat.location.parts[-1].end.position  #Featureの終了位置
       	strand = feat.strand
        if feat.type == "CDS":
            if "gene" in feat.qualifiers:
                gene = feat.qualifiers["gene"][0]
            else:
                gene = "N.A."
            
            if "product" in feat.qualifiers:
                product = feat.qualifiers["product"][0]
                if "," in product:
                    product = product.replace(",",":")
            else:
                product = "N.A."
            
            if "protein_id" in feat.qualifiers:
                    protein_id = feat.qualifiers["protein_id"][0]
            else:
                protein_id = "N.A."
            results.append([gene,protein_id,str(start),str(end),str(strand),product])  
    return results

if __name__ == "__main__":
    for gbk in SeqIO.parse("../sample_genome/GCF_000005845.2_ASM584v2_genomic.gbff","genbank"):
        results = get_CDS_info(gbk)
        with open(gbk.id + "_CDS_info.csv","w") as o:
            for result in results:
                o.write(",".join(result) + "\n") 

