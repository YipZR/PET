from Bio import Entrez
from Bio import SeqIO
import json
import xmltodict
import os


os.chdir('D:/PET/invertebrate')

Entrez.email = "122616990@qq.com"
handle = Entrez.esearch(db="sra", term="invertebrate", retmax="334819")
record = Entrez.read(handle)
handle.close()
id_list = record["IdList"]

for i in range(0,len(id_list)):
    try:
        handle_2 = Entrez.efetch(db="sra",id=id_list[i])
        xml = handle_2.read()
        json_file="invertebrate.json"
        convertJson = xmltodict.parse(xml,encoding = 'utf-8')
        jsonStr = json.dumps(convertJson,indent=1)

        with open(json_file, 'w+',encoding = 'utf-8') as f:
            f.write(jsonStr)

        with open('invertebrate.json','r') as f:
            data = json.load(f)
        file=open('invertebrate.csv','a')
        
        Run_acc = data['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['RUN_SET']['RUN']['@accession']
        Bioproject = data['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['STUDY']['IDENTIFIERS']['EXTERNAL_ID']['#text']
        Biosamples = data['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['SAMPLE']['IDENTIFIERS']['EXTERNAL_ID']['#text']
        Organism = data['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['SAMPLE']['SAMPLE_NAME']['SCIENTIFIC_NAME']
        Strategy = data['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['EXPERIMENT']['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_STRATEGY']
        Source =  data['EXPERIMENT_PACKAGE_SET']['EXPERIMENT_PACKAGE']['EXPERIMENT']['DESIGN']['LIBRARY_DESCRIPTOR']['LIBRARY_SOURCE']

        item=','.join([Run_acc,Bioproject,Biosamples,Strategy,Organism,Source])
        file.writelines(item+'\n')
        file.closes
    except:
        continue
