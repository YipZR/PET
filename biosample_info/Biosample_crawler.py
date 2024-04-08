from urllib import request
import csv
from bs4 import BeautifulSoup

#https://www.ncbi.nlm.nih.gov/biosample/SAMN27259107

Env='MHETuniq'
file = open(f'/home/yzr/Work/Meta/biosample_info/{Env}_biosample.txt')
for line in file:
    Biosample = line.strip()
    url=f'https://www.ncbi.nlm.nih.gov/biosample/{Biosample}'
    req=request.Request(url,headers={'User-Agent':'Mozilla/5.0 (compatible; MSIE 5.5; Windows NT)'})
    page=request.urlopen(req)
    bf=BeautifulSoup(page, "html.parser")
    page.close()

    output = open(f"{Env}_parser.txt","a")
    
    try:
        output.write(Biosample + '\t')
        table=bf.findAll('table',{'class':'docsum'})
        tab=table[0]
        for tr in tab.find_all('tr'):
            for td in tr.find_all('td'):
                output.write(td.text + '\t')
        output.write('\n')
        output.close()
    except:
        continue
file.close()
