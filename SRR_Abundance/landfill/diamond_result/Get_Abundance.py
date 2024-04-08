from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio import Entrez
import pandas as pd
import numpy as np
import os,re
def get_genes_abundance(file,all_gene,sra_id):
    diamond_data=pd.read_table(path+file,sep='\t',encoding='utf-8',quoting=3,header=None)
    diamond_data.columns=["query","hit","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore"]
    diamond_data=diamond_data.loc[diamond_data['pident']>=80,:]
    diamond_data=pd.merge(diamond_data,all_gene,how="left",left_on="hit",right_on="description")
    gene_count=pd.DataFrame(diamond_data["gene_name"].value_counts(sort=False))
    gene_count.columns=["count"]
    gene_count["gene_name"] = gene_count.index
    gene_count=gene_count.reset_index(drop=True)

    gene_length=pd.DataFrame(diamond_data.groupby("gene_name")["kb"].median())
    gene_length["gene_name"] = gene_length.index
    gene_length=gene_length.reset_index(drop=True)

    gene_count=pd.merge(gene_count,gene_length,how="left",on="gene_name")
    gene_count["rpk"]=gene_count[["count","kb"]].apply(lambda x:x["count"]/x["kb"],axis=1)
    all_gene_categories=pd.read_table("/home/yzr/Work/Meta/MetaQ/categories_merged.txt",sep='\t',encoding='utf-8',quoting=3,header=0)  
    gene_count=pd.merge(all_gene_categories,gene_count,how="left",on="gene_name")
    GC1=gene_count[gene_count["category"]=="reference"]
    GC1=pd.DataFrame(GC1.groupby("category")["rpk"].median())
    GC2=gene_count[gene_count["category"]!="reference"]
    GC2=pd.DataFrame(GC2.groupby("category")["rpk"].sum())
    gene_cluster_count= pd.concat([GC1,GC2])
    #拆成两部分，ref和非ref，ref做median，非ref做sum，然后合并，后面照旧
    #gene_cluster_count=pd.DataFrame(gene_count.groupby("category")["rpk"].median())
    gene_cluster_count.columns=["count"]
    gene_cluster_count["category"] = gene_cluster_count.index
    gene_cluster_count=gene_cluster_count.reset_index(drop=True)
    #print(gene_cluster_count)
    #gene_cluster_count.to_csv("result_.txt",sep="\t")

    ref=gene_cluster_count[gene_cluster_count["category"]=="reference"]
    ref=float(ref["count"])
    gene_cluster_count[sra_id]=gene_cluster_count["count"]/ref*1000
    #print(gene_cluster_count)
    del gene_cluster_count["count"]
    return (gene_cluster_count)


def get_PGPR_database(seq):
    IDs=[]
    sequence=[]
    for seq_record in SeqIO.parse("/home/yzr/Work/Meta/MetaQ/merged_MetaQ.faa","fasta"):
        IDs.append(str(seq_record.description))
        sequence.append(str(seq_record.seq))
    all_gene=pd.DataFrame(data={"description":IDs,"sequence":sequence})
    all_gene['kb']=all_gene['sequence'].str.len()/1000
    df=all_gene['description'].str.split('|',expand=True)    
    df.columns=["gene_name","significance","metabolites","gene_id","genome_id","saa"]
    all_gene=all_gene.join(df["gene_name"])
    return (all_gene)
seq="/home/yzr/Work/Meta/MetaQ/merged_MetaQ.faa"
#all_gene.to_csv("result_.txt",sep="\t")
all_gene=get_PGPR_database(seq)
path="/data_1/yzr_data/SRR/landfill/diamond_result/"
files=os.listdir(path)
os.chdir("/data_1/yzr_data/SRR/landfill/diamond_result/")

result=pd.read_table("/home/yzr/Work/Meta/MetaQ/categories_merged.txt",sep='\t',encoding='utf-8',quoting=3,header=0)
result=pd.DataFrame(result["category"])
result=result.drop_duplicates()
#print(result)

for file in files:
    if re.search('.m8',file):
        ID=file.split("_")[0]
        sra_id=ID
        print(sra_id)
        gene_adundance=get_genes_abundance(file,all_gene,sra_id)
        result=pd.merge(result,gene_adundance,how="left",on="category")

result=pd.DataFrame(result.values.T,index=result.columns,columns=result.index)
array=np.array(result)
list_=array.tolist()
list_=list_[0]
result.columns=list_
result.drop(result.index[0],inplace=True)

sra_sample_table=pd.read_table("/home/yzr/Work/Meta/EBI/landfill_ID_SRA.txt",sep='\t',encoding='utf-8',quoting=3,header=0)

result=result.join(sra_sample_table)

result.to_csv("result.txt",sep="\t")


