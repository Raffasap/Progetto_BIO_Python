#!/usr/bin/env python
# coding: utf-8

# In[1]:


import re


# In[2]:


gtf_file_name = './input.gtf'


# In[3]:


with open(gtf_file_name, 'r') as gtf_input_file:
    gtf_file_rows = gtf_input_file.readlines()


# In[4]:


# estrazione dei primi due elementi della prima riga per controllare 
# successivamente che tutti i record contengano lo stesso "seqname" e la stessa "source" 

nome = gtf_file_rows[0].split('\t')[0]
sorgente = gtf_file_rows[0].split('\t')[1]


# In[5]:


#definisco l'array che andrà a contenere le violazioni presenti nel file .gtf e le funzioni di verifica dei signoli valori di ogni record
violations = []

def validate_num_fields(row, index): # il numero di elementi di ogni record non può essere minore di 9 o maggiore di 10 (considerando il campo "commenti" opzionale)
    fields = len(row.split('\t'))
    if fields < 9 or fields > 10:
        violations.append('record: ' + str(index) + ' - ogni record deve contenere almeno 9 e non più di 10 elementi, ognuno separato dal carattere tab')
        return False
    else: 
        return True

def validate_seqname(value, index): # tutti i record dello stesso file devono avere "seqname" univoco
    if value != nome : 
        return violations.append('record: ' + str(index) + ' - il campo "seqname" deve essere uguale per tutti i record')

def validate_source(value, index): # tutti i record dello stesso file devono avere "source" univoco
    if value != sorgente : 
        return violations.append('record: ' + str(index) + ' - il campo "source" deve essere uguale per tutti i record')

def validate_feature(value, index): # valida l'elemento "feature" e conteggio del numero di record "required"    
    global CDS_count, start_codon_count, stop_codon_count, exon_count
    if value not in ['exon', 'CDS', '5UTR', '3UTR', 'start_codon', 'stop_codon', 'inter', 'inter_CNS', 'intron_CNS']:
        return violations.append('record: ' + str(index) + ' - valore del campo "feature" non ammesso')
    elif value == 'CDS':
        CDS_count += 1
    elif value == 'start_codon':
        start_codon_count += 1
    elif value == 'stop_codon':
        stop_codon_count += 1
    elif value == 'exon':
        exon_count += 1

def validate_start(value, index): # l'elemento "start" deve avere come valore un numero intero > 0
    if value == '0':
        return violations.append('record: ' + str(index) + ' - il valore del campo "start" deve essere un valore numerico intero non negativo')
    elif value.isnumeric() == False:
        return violations.append('record: ' + str(index) + ' - il valore del campo "start" deve essere un valore numerico intero non negativo')

def validate_end(value, index): # l'elemento "end" deve avere come valore un numero intero > 0, "end" non può essere minore di "start" e inoltre se la "feature" è "start_codon"/"stop_codon" la lunghezza massima è 3bp
    if value == '0':
         return violations.append('record: ' + str(index) + ' - il valore del campo "end" deve essere un valore numerico intero non negativo')
    elif value.isnumeric() == False:
         return violations.append('record: ' + str(index) + ' - il valore del campo "end" deve essere un valore numerico intero non negativo')
    elif int(end) < int(start):
         return violations.append('record: ' + str(index) + ' - il valore del campo "end" deve essere maggiore o uguale al campo "start"')
    elif feature in ['stop_codon', 'start_codon'] and (int(end) - int(start)) > 3:
         return violations.append('record: ' + str(index) + ' - le funzioni "start_codon" e "stop_codon" devono essere lunghe al massimo 3 bp')

def validate_score(value, index): # se il valore è diverso da '.' deve essere un intero o un float, dato che int è sottoinsieme di float basta controllare che sia float
    if value != '.':
        try:
            float(value)
        except:
            violations.append('record: ' + str(index) + ' - se il valore del campo "score" è diverso da "." deve essere un intero o un numero in virgola mobile')

def validate_strand(value, index): # consentiti solo i valori +/-
    if value not in ['-', '+']:
        return violations.append('record: ' + str(index) + ' - valore del campo "strand" diverso da quelli consentiti')
    return

def validate_frame(value, index): # il campo "frame" può assumere valore (0, 1, 2) solo se "feature" è tra (CDS, start_codon, stop_codon) altrimenti deve valere "."
    if feature in ['CDS', 'start_codon', 'stop_codon']:    
        if value not in ['0', '1', '2']:
            return violations.append('record: ' + str(index) + ' - valore del campo "frame" diverso da quelli consentiti(0, 1, 2)')
    else:
        if frame != '.':
            return violations.append('record: ' + str(index) + ' - il campo frame deve assumere valore "." se non si tratta di una delle seguenti features: CDS, start_codon, stop_codon')
    
def validate_attributes(value, index): # l'elemento "attributes" deve contenere per forza il transcript_id e il gene_id, gli altri attributi non vengono calcolati
    gene_id = re.findall('gene_id\s+"(\w+)";', value)
    transcript_id = re.findall('transcript_id\s+"([^"]+)";', value)
    if not transcript_id :   
        return violations.append('record: ' + str(index) + ' - il campo "attributes" deve contenere il "transcript_id"')
    if not gene_id :   
        return violations.append('record: ' + str(index) + ' - il campo "attributes" deve contenere il "gene_id"')

def validate_count(CDS, start, stop, exon): # "feature" required
    if CDS == 0:
        violations.append('il file deve contenere almeno un record con feature "CDS"')
    if start == 0:
        violations.append('il file deve contenere almeno un record con feature "start_codon"')
    if stop == 0:
        violations.append('il file deve contenere almeno un record con feature "stop_codon"')
    if exon == 0:
        violations.append('il file deve contenere almeno un record con feature "exon"')


# In[6]:


i = 0
CDS_count = 0
start_codon_count = 0
stop_codon_count = 0
exon_count = 0

for row in gtf_file_rows: # ciclo sulla lista che contiene i record del file, splitto ogni record nei suoi elementi, verifico che il numero di elementi sia giusto, e per ognuno di essi invoco la propria funzione validante
    
    if validate_num_fields(row, i): 
       
        seqname = row.split('\t')[0]
        source = row.split('\t')[1]
        feature = row.split('\t')[2]
        start = row.split('\t')[3]
        end = row.split('\t')[4]
        score = row.split('\t')[5]
        strand = row.split('\t')[6]
        frame = row.split('\t')[7]
        attributes = row.split('\t')[8]
        
        validate_seqname(seqname, i)
        validate_source(source, i)
        validate_feature(feature, i)
        validate_start(start, i)
        validate_end(end, i)
        validate_score(score, i)
        validate_strand(strand, i)
        validate_frame(frame, i)
        validate_attributes(attributes, i)

    i = i + 1

validate_count(CDS_count, start_codon_count, stop_codon_count, exon_count) # dopo aver controllato tutti i record verifico la presenza delle "feature" required


# In[7]:


for violation in violations: # produco in output l'elenco delle violazioni presenti nel file
    print(violation) 

