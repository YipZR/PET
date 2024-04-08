from Bio import Entrez
from Bio import SeqIO
import json
import xmltodict
import os


os.chdir('D:/PET/biosample')

Entrez.email = "122616990@qq.com"
id_list = 'SAMN33826073'

for i in range(0,len(id_list)):
    try:
        handle = Entrez.efetch(db="BioSample",id='SAMN33826073')
        xml = handle.read()
        json_file="biosample.json"
        convertJson = xmltodict.parse(xml,encoding = 'utf-8')
        jsonStr = json.dumps(convertJson,indent=1)

        with open(json_file, 'w+',encoding = 'utf-8') as f:
            f.write(jsonStr)

        with open('biosample.json','r') as f:
            data = json.load(f)
        file=open('biosample.csv','a')
        
        Biosamples = data['BioSampleSet']['BioSample']['@accession']
        Run_acc = data['BioSampleSet']['BioSample']['Ids']['Id']['2']['#text']
        Description = data['BioSampleSet']['BioSample']['Description']['Comment']['Paragraph']
        Source =  data['BioSampleSet']['BioSample']['Attributes']['Attribute']

        item=','.join([Run_acc,Biosamples,Description])
        file.writelines(item+'\n')
        file.closes
    except:
        continue
