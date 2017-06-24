import sys
from Bio import Entrez

email = 'fe.sepulveda.c@gmail.com'

Entrez.email = email

# sample = sys.argv[1]


# handle = Entrez.esearch('BioSample', term=sample, id='text')
# record = Entrez.read(handle)

# #print record

# id_sample = record['IdList'][0] 
# handle = Entrez.efetch('BioSample', id=id_sample, retmode='text')



# date = ''
# disease = ''
# geographic = ''
# isolation = ''


# for line in handle.readlines():
# 	if 'collection date' in line:
# 		date = line.split('"')[1]
# 		#print date
# 	if 'host disease' in line:
# 		disease = line.split('"')[1]
# 		#print disease
# 	if 'geographic location' in line:
# 		geographic = line.split('"')[1]
# 		#print geographic
# 	if 'isolation source' in line:
# 		isolation = line.split('"')[1]
# 		#print isolation

# if not date:
# 	date = 'NN'
# if not disease:
# 	disease = 'NN'
# if not geographic:
# 	geographic = 'NN'
# if not isolation:
# 	isolation = 'NN'

# print(sample, date, disease, isolation, geographic)
from Bio import SeqIO
sample = sys.argv[1]
handle = Entrez.efetch(db='nuccore', id=sample, rettype="gb", retmode='text')

record = SeqIO.read(handle, 'genbank')

print sample, record.dbxrefs