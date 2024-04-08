import pandas as pd
import warnings
warnings.filterwarnings("ignore")
#Kingdom Phylum Class Order Family Genus Species
df = pd.read_table("Run_OTU_result.txt")
#print(df1)

clade_split = df['clade_name'].str.split('|')

Kingdom = df[clade_split.apply(lambda x: x[-1].startswith('d__'))]
Kingdom_to_normalize = Kingdom.columns.difference(['clade_name'])
Kingdom[Kingdom_to_normalize] = Kingdom[Kingdom_to_normalize].apply(lambda x: x / x.sum() * 100, axis=0)
Kingdom.to_csv("OTU_Kingdom.txt",sep="\t",index=False)

Phylum = df[clade_split.apply(lambda x: x[-1].startswith('p__'))]
Phylum_to_normalize = Phylum.columns.difference(['clade_name'])
Phylum[Phylum_to_normalize] = Phylum[Phylum_to_normalize].apply(lambda x: x / x.sum() * 100, axis=0)
Phylum.to_csv("OTU_Phylum.txt",sep="\t",index=False)

Order = df[clade_split.apply(lambda x: x[-1].startswith('o__'))]
Order_to_normalize = Order.columns.difference(['clade_name'])
Order[Order_to_normalize] = Order[Order_to_normalize].apply(lambda x: x / x.sum() * 100, axis=0)
Order.to_csv("OTU_Order.txt",sep="\t",index=False)

Family = df[clade_split.apply(lambda x: x[-1].startswith('f__'))]
Family_to_normalize = Family.columns.difference(['clade_name'])
Family[Family_to_normalize] = Family[Family_to_normalize].apply(lambda x: x / x.sum() * 100, axis=0)
Family.to_csv("OTU_Family.txt",sep="\t",index=False)

Genus = df[clade_split.apply(lambda x: x[-1].startswith('g__'))]
Genus_to_normalize = Genus.columns.difference(['clade_name'])
Genus[Genus_to_normalize] = Genus[Genus_to_normalize].apply(lambda x: x / x.sum() * 100, axis=0)
Genus.to_csv("OTU_Genus.txt",sep="\t",index=False)

Species = df[clade_split.apply(lambda x: x[-1].startswith('s__'))]
Species_to_normalize = Species.columns.difference(['clade_name'])
Species[Species_to_normalize] = Species[Species_to_normalize].apply(lambda x: x / x.sum() * 100, axis=0)
Species.to_csv("OTU_Species.txt",sep="\t",index=False)
