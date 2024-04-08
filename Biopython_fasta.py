import pandas as pd
from Bio import Entrez,SeqIO,Seq
Entrez.email = 'yezr20@lzu.edu.cn'
ID_list = pd.read_table('D:/PET_Fasta/PETorthologs/qcover80ID.txt')
list = []
for i in range(0,len(ID_list)):
    with open('D:/PET_Fasta/PETorthologs/qcover80ID.txt') as f:
        line = f.readlines()[i].rstrip()
        handle_fa = Entrez.efetch(db="protein", id=line, rettype="fasta", retmode="text")
        read_fa = handle_fa.read()
        with open('D:/PET_Fasta/PETorthologs/qcover80.fasta',"a")as file:
            file.write(read_fa)