import argparse
#解析命令行参数和选项的标准模块
import os
import sys
import re
import pandas as pd
from itertools import takewhile


def merge(file_list, out_file):

    merged_tables = pd.DataFrame()
    for file in file_list:
        iIn = pd.read_csv(file,
                          sep='\t',
                          names=["clade_name","relative_abundance"],
                          index_col=0
                          )
        #print(iIn)
        if merged_tables.empty:
            merged_tables = iIn.iloc[:, 0].rename(os.path.splitext(os.path.basename(file))[0]).to_frame()
        else:
            merged_tables = pd.merge(iIn.iloc[:, 0].rename(os.path.splitext(os.path.basename(file))[0]).to_frame(),
                                     merged_tables,
                                     how='outer',
                                     left_index=True,
                                     right_index=True
                                     )

    merged_tables.fillna('0').reset_index().to_csv(out_file, index=False, sep='\t')


path="/data_1/yzr_data/SRR/Not/kraken2"
files=os.listdir(path)
os.chdir(path)

file_list=[]
for file in files:
    if re.search('_kraken2.txt',file):
        file_list.append(file)

with open("Not_Abundance_OTU_result.txt", 'w') as out_file:
    merge(file_list, out_file)
