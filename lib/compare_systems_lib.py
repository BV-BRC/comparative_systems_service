#!/usr/bin/env python

import copy
import gzip
import json
import multiprocessing
import os
import re
import shutil
import subprocess
import sys
import tarfile
import urllib.request as request
from contextlib import closing
from multiprocessing import Process

import requests
import pandas as pd
import numpy as np

from bvbrc_api import authenticateByEnv,getGenomeIdsByGenomeGroup,getFeatureDataFrame,getSubsystemsDataFrame,getPathwayDataFrame,getDataForGenomes,getQueryData,getQueryDataText

import time
import io

def chunker(seq, size):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def build_key(component_list):
    mod_list = []
    for x in component_list:
        if not isinstance(x,str):
            x = 'nan'
        mod_list.append(x)
    #key_parts = [x.lower() for x in mod_list]
    return ":".join(mod_list)

# Returns the entry with the maximum occurences in a pandas dataframe column
# df: pandas dataframe
# col: column name
def get_maximum_value(df, col):
    # value counts returns an ordered list, descending
    col_counts = df[col].value_counts() 
    max_label = col_counts[0].index
    return max_label

def return_columns_to_remove(system,columns):
    if system is 'subsystems_genes':
        drop_subsystem_columns = ['date_inserted','date_modified','genome_name','gene','owner','patric_id','public','product','refseq_locus_tag','taxon_id','_version_']
        table_columns = list(set.intersection(set(drop_subsystem_columns),set(columns))) 
        return table_columns
    elif system is 'subsystems_subsystems':
        drop_subsystem_columns = ['feature_id','public','role_id','genome_id','taxon_id','role_name','owner','product','patric_id','genome_name','id','_version_','date_inserted','date_modified']
        table_columns = list(set.intersection(set(drop_subsystem_columns),set(columns))) 
        return table_columns
    elif system is 'proteinfamilies_plfams':
        drop_plfams_columns = ['genome_name','accession','patric_id','refseq_locus_tag',
                                    'alt_locus_tag','feature_id','annotation','feature_type',
                                    'start','end','strand','figfam_id','pgfam_id','protein_id',
                                    'aa_length','gene','go']
        table_columns = list(set.intersection(set(drop_plfams_columns),set(columns))) 
        return table_columns
    elif system is 'proteinfamilies_pgfams':
        drop_pgfams_columns = ['genome_name','accession','patric_id','refseq_locus_tag',
                                    'alt_locus_tag','feature_id','annotation','feature_type',
                                    'start','end','strand','figfam_id','plfam_id','protein_id',
                                    'aa_length','gene','go']
        table_columns = list(set.intersection(set(drop_pgfams_columns),set(columns))) 
        return table_columns
    elif system is 'pathways_genes':
        drop_gene_columns = ['genome_name','accession','alt_locus_tag','refseq_locus_tag','feature_id','annotation','product']
        table_columns = list(set.intersection(set(drop_gene_columns),set(columns)))
        return table_columns
    else: # pathways does not have drop columns
        sys.stderr.write("Error, system is not a valid type\n")
        return [] 

def run_families_v2(genome_ids, query_dict, output_file, output_dir, genome_data, session):
    data_dict = {} 
    data_dict['plfam'] = {}
    data_dict['pgfam'] = {}
    plfam_genomes = {}
    pgfam_genomes = {}
    present_genome_ids = set()
    genomes_missing_data = {}
    for gids in chunker(genome_ids, 20):
        base = "https://www.patricbrc.org/api/genome_feature/?http_download=true"
        query = f"in(genome_id,({','.join(gids)}))&limit(2500000)&sort(+feature_id)&eq(annotation,PATRIC)"
        headers = {"accept":"text/tsv", "content-type":"application/rqlquery+x-www-form-urlencoded", 'Authorization': session.headers['Authorization']}
        result_header = True
        for line in getQueryData(base,query,headers):
            if result_header:
                result_header = False
                print(line)
                continue
            line = line.strip().split('\t')
            # 20 entries in query result with pgfam and plfam data
            if len(line) < 20: 
                continue
            try:
                genome_id = line[1].replace('\"','')
                plfam_id = line[14].replace('\"','')
                pgfam_id = line[15].replace('\"','')
                aa_length = line[17].replace('\"','')
                product = line[19].replace('\"','')
            except Exception as e:
                sys.stderr.write(f'Error with the following line:\n{e}\n{line}\n')
                continue
            if aa_length == '':
                continue
            present_genome_ids.add(genome_id)
            ### add to missing genomes data dict
            if genome_id not in genomes_missing_data:
                genomes_missing_data[genome_id] = True
            ### plfam counts
            if plfam_id not in data_dict['plfam']:
                data_dict['plfam'][plfam_id] = {} 
                data_dict['plfam'][plfam_id]['aa_length_list'] = [] 
                data_dict['plfam'][plfam_id]['feature_count'] = 0 
                data_dict['plfam'][plfam_id]['genome_count'] = 0 
                data_dict['plfam'][plfam_id]['product'] = product 
            if plfam_id not in plfam_genomes:
                plfam_genomes[plfam_id] = {} 
            if genome_id not in plfam_genomes[plfam_id]:
                plfam_genomes[plfam_id][genome_id] = 0
            data_dict['plfam'][plfam_id]['aa_length_list'].append(int(aa_length))
            data_dict['plfam'][plfam_id]['feature_count']+=1
            data_dict['plfam'][plfam_id]['genome_count'] = len(plfam_genomes[plfam_id])
            plfam_genomes[plfam_id][genome_id]+=1
            ### pgfam counts
            if pgfam_id not in data_dict['pgfam']:
                data_dict['pgfam'][pgfam_id] = {} 
                data_dict['pgfam'][pgfam_id]['aa_length_list'] = [] 
                data_dict['pgfam'][pgfam_id]['feature_count'] = 0 
                data_dict['pgfam'][pgfam_id]['genome_count'] = 0 
                data_dict['pgfam'][pgfam_id]['product'] = product 
            if pgfam_id not in pgfam_genomes:
                pgfam_genomes[pgfam_id] = {} 
            if genome_id not in pgfam_genomes[pgfam_id]:
                pgfam_genomes[pgfam_id][genome_id] = 0
            data_dict['pgfam'][pgfam_id]['aa_length_list'].append(int(aa_length))
            data_dict['pgfam'][pgfam_id]['feature_count']+=1
            data_dict['pgfam'][pgfam_id]['genome_count'] = len(pgfam_genomes[pgfam_id])
            pgfam_genomes[pgfam_id][genome_id]+=1

    # TODO: https://alpha.bv-brc.org/api/protein_family_ref/
    # - get protein family description data
    product_dict = {}
    for plids_list in chunker(list(data_dict['plfam'].keys()),5000):
        base = "https://alpha.bv-brc.org/api/protein_family_ref/?http_download=true"
        query = f"in(family_id,({','.join(plids_list)}))&limit(2500000)&sort(+family_id)"
        headers = {"accept":"application/json", "content-type":"application/rqlquery+x-www-form-urlencoded", 'Authorization': session.headers['Authorization']}
        #headers = {"accept":"text/tsv", "content-type":"application/rqlquery+x-www-form-urlencoded", 'Authorization': session.headers['Authorization']}
        res_data = getQueryDataText(base,query,headers)
        text_data = json.loads(res_data)
        for entry in text_data:
            product_dict[entry['family_id']] = entry['family_product']
    for pgids_list in chunker(list(data_dict['pgfam'].keys()),5000):
        base = "https://alpha.bv-brc.org/api/protein_family_ref/?http_download=true"
        query = f"in(family_id,({','.join(pgids_list)}))&limit(2500000)&sort(+family_id)"
        headers = {"accept":"application/json", "content-type":"application/rqlquery+x-www-form-urlencoded", 'Authorization': session.headers['Authorization']}
        res_data = getQueryDataText(base,query,headers)
        text_data = json.loads(res_data)
        for entry in text_data:
            product_dict[entry['family_id']] = entry['family_product']

    # go back and get the mean, max, min, std dev for each family_id
    plfam_line_list = []        
    pgfam_line_list = []
    plfam_genome_list = {}
    pgfam_genome_list = {}
    genome_str_dict = {}
    genome_str_dict['plfam'] = {}
    genome_str_dict['pgfam'] = {}
    for plfam_id in data_dict['plfam']:
        if plfam_id == '':
            continue
        aa_length_list = data_dict['plfam'][plfam_id]['aa_length_list']
        aa_length_max = max(aa_length_list)
        aa_length_min = min(aa_length_list)
        aa_length_mean = np.mean(aa_length_list)
        aa_length_std = np.std(aa_length_list)
        feature_count = data_dict['plfam'][plfam_id]['feature_count']
        genome_count = data_dict['plfam'][plfam_id]['genome_count']
        #genomes = format(feature_count,'#04x').replace('0x','')
        genomes_dir = {} 
        plfam_genome_list[plfam_id] = []
        for gid in genome_ids:
            if gid in plfam_genomes[plfam_id]:
                genomes_dir[gid] = format(plfam_genomes[plfam_id][gid],'#04x').replace('0x','')
                plfam_genome_list[plfam_id].append(gid)
            else:
                genomes_dir[gid] = '00'
        #genomes = ''.join(genomes_list)
        genome_str_dict['plfam'][plfam_id] = genomes_dir 
        #product = data_dict['plfam'][plfam_id]['product']
        if plfam_id in product_dict:
            product = product_dict[plfam_id]
        else:
            print('{0} not in product_dict'.format(plfam_id))
            product = 'NOTHING'
        plfam_str = f'{plfam_id}\t{feature_count}\t{genome_count}\t{product}\t{aa_length_min}\t{aa_length_max}\t{aa_length_mean}\t{aa_length_std}'
        plfam_line_list.append(plfam_str)
    for pgfam_id in data_dict['pgfam']:
        if pgfam_id == '':
            continue
        aa_length_list = data_dict['pgfam'][pgfam_id]['aa_length_list']
        aa_length_max = max(aa_length_list)
        aa_length_min = min(aa_length_list)
        aa_length_mean = np.mean(aa_length_list)
        aa_length_std = np.std(aa_length_list)
        feature_count = data_dict['pgfam'][pgfam_id]['feature_count']
        genome_count = data_dict['pgfam'][pgfam_id]['genome_count']
        #genomes = format(feature_count,'#04x').replace('0x','')
        genomes_dir = {}
        pgfam_genome_list[pgfam_id] = []
        for gid in genome_ids:
            if gid in pgfam_genomes[pgfam_id]:
                #genomes+=format(pgfam_genomes[pgfam_id][gid],'#04x').replace('0x','')
                genomes_dir[gid] = format(pgfam_genomes[pgfam_id][gid],'#04x').replace('0x','')
                pgfam_genome_list[pgfam_id].append(gid)
            else:
                genomes_dir[gid] = '00'
        #genomes = ''.join(genomes_list)
        genome_str_dict['pgfam'][pgfam_id] = genomes_dir 
        #product = data_dict['pgfam'][pgfam_id]['product']
        product = product_dict[pgfam_id]
        pgfam_str = f'{pgfam_id}\t{feature_count}\t{genome_count}\t{product}\t{aa_length_min}\t{aa_length_max}\t{aa_length_mean}\t{aa_length_std}'
        pgfam_line_list.append(pgfam_str)

    #output_json['genome_ids'] = genome_ids
    #output_json['genome_ids'] = list(set(genome_ids).intersection(present_genome_ids)) 

    unsorted_genome_ids = [gid for gid in genome_ids if gid in present_genome_ids] 
    tmp_data = genome_data.loc[genome_data['Genome ID'].isin(unsorted_genome_ids)]
    tmp_data.set_index('Genome ID',inplace=True)
    tmp_data = tmp_data.loc[unsorted_genome_ids]
    unsorted_genome_names = tmp_data['Genome Name'].tolist()
    sorted_genome_names, sorted_genome_ids = zip(*sorted(zip(unsorted_genome_names,unsorted_genome_ids)))

    # add genomes string to each line
    for x in range(0,len(plfam_line_list)): 
        line_parts = plfam_line_list[x].split('\t')
        genomes_dir = genome_str_dict['plfam'][line_parts[0]]
        genome_str = ''
        for gid in sorted_genome_ids:
            genome_str+=genomes_dir[gid]
        line_parts.append(genome_str)
        plfam_line_list[x] = '\t'.join(line_parts)
    for x in range(0,len(pgfam_line_list)): 
        line_parts = pgfam_line_list[x].split('\t')
        genomes_dir = genome_str_dict['pgfam'][line_parts[0]]
        genome_str = ''
        for gid in sorted_genome_ids:
            genome_str+=genomes_dir[gid]
        line_parts.append(genome_str)
        pgfam_line_list[x] = '\t'.join(line_parts)

    
    header = 'family_id\tfeature_count\tgenome_count\tproduct\taa_length_min\taa_length_max\taa_length_mean\taa_length_std\tgenomes'
    plfam_line_list.insert(0,header)
    pgfam_line_list.insert(0,header)

    output_json = {}
    output_json['plfam'] = '\n'.join(plfam_line_list) 
    output_json['pgfam'] = '\n'.join(pgfam_line_list) 
    output_json['genome_ids'] = sorted_genome_ids 
    output_json['genome_names'] = sorted_genome_names
    output_json['job_name'] = output_file
    output_json['plfam_genomes'] = plfam_genome_list 
    output_json['pgfam_genomes'] = pgfam_genome_list 

    output_json_file = os.path.join(output_dir,output_file+'_proteinfams_tables.json')
    with open(output_json_file,"w") as o:
        o.write(json.dumps(output_json))

    print("ProteinFamilies Complete")
    return ({ 'success': True, 'genomes': present_genome_ids })

def run_families(genome_ids, query_dict, output_file, output_dir, genome_data, session):
    plfam_dict = {}
    pgfam_dict = {}
    genome_dict = {}
    plfam_dict['unique_set'] = set()
    pgfam_dict['unique_set'] = set()
    present_genome_ids = set()
    for gids in chunker(genome_ids, 20):
        base = "https://www.patricbrc.org/api/genome_feature/?http_download=true"
        query = f"in(genome_id,({','.join(gids)}))&limit(2500000)&sort(+feature_id)&eq(annotation,PATRIC)"
        headers = {"accept":"text/tsv", "content-type":"application/rqlquery+x-www-form-urlencoded", 'Authorization': session.headers['Authorization']}
        result_header = True
        for line in getQueryData(base,query,headers):
            if result_header:
                result_header = False
                print(line)
                continue
            line = line.strip().split('\t')
            # 20 entries in query result with pgfam and plfam data
            if len(line) < 20: 
                continue
            try:
                genome_id = line[1].replace('\"','')
                plfam_id = line[14].replace('\"','')
                pgfam_id = line[15].replace('\"','')
                aa_length = line[17].replace('\"','')
                product = line[19].replace('\"','')
            except Exception as e:
                sys.stderr.write(f'Error with the following line:\n{e}\n{line}\n')
                continue
            if genome_id not in genome_dict:
                genome_dict[genome_id] = {}
                genome_dict[genome_id]['plfam_set'] = set()
                genome_dict[genome_id]['pgfam_set'] = set()
            present_genome_ids.add(genome_id)
            genome_dict[genome_id]['plfam_set'].add(plfam_id) 
            genome_dict[genome_id]['pgfam_set'].add(pgfam_id) 
            if genome_id not in plfam_dict:
                plfam_dict[genome_id] = {}
            if genome_id not in pgfam_dict:
                pgfam_dict[genome_id] = {}
            if plfam_id and plfam_id not in plfam_dict[genome_id]:
                plfam_dict[genome_id][plfam_id] = {}
                plfam_dict[genome_id][plfam_id]['aa_length_list'] = []
                # TODO: check if I need to check for duplicate features
                plfam_dict[genome_id][plfam_id]['feature_count'] = 0
                plfam_dict[genome_id][plfam_id]['genome_count'] = 1
                plfam_dict[genome_id][plfam_id]['product'] = product
            if pgfam_id and pgfam_id not in pgfam_dict[genome_id]:
                pgfam_dict[genome_id][pgfam_id] = {}
                pgfam_dict[genome_id][pgfam_id]['aa_length_list'] = []
                # TODO: check if I need to check for duplicate features
                pgfam_dict[genome_id][pgfam_id]['feature_count'] = 0
                pgfam_dict[genome_id][pgfam_id]['genome_count'] = 1
                pgfam_dict[genome_id][pgfam_id]['product'] = product
            if plfam_id:
                plfam_dict['unique_set'].add(plfam_id)
                plfam_dict[genome_id][plfam_id]['aa_length_list'].append(int(aa_length))
                plfam_dict[genome_id][plfam_id]['feature_count'] = plfam_dict[genome_id][plfam_id]['feature_count'] + 1
            if pgfam_id:
                pgfam_dict['unique_set'].add(pgfam_id)
                pgfam_dict[genome_id][pgfam_id]['aa_length_list'].append(int(aa_length))
                pgfam_dict[genome_id][pgfam_id]['feature_count'] = pgfam_dict[genome_id][pgfam_id]['feature_count'] + 1
        
    plfam_line_list = []        
    pgfam_line_list = []
    header = 'family_id\tgenome_id\tfeature_count\tgenome_count\tproduct\taa_length_min\taa_length_max\taa_length_mean\taa_length_std\tgenomes'
    plfam_line_list.append(header)
    pgfam_line_list.append(header)
    for plfam_id in plfam_dict['unique_set']: 
        for gid in genome_ids: 
            if gid not in plfam_dict:
                continue
            #plfam_data['plfam_id'] = plfam_id
            #plfam_data['genome_id'] = gid
            if plfam_id in plfam_dict[gid]:
                aa_length_list = plfam_dict[gid][plfam_id]['aa_length_list'] 
                aa_length_max = max(aa_length_list)
                aa_length_min = min(aa_length_list)
                aa_length_mean = np.mean(aa_length_list)
                aa_length_std = np.std(aa_length_list)
                feature_count = plfam_dict[gid][plfam_id]['feature_count']
                genome_count = plfam_dict[gid][plfam_id]['genome_count']
                genomes = format(feature_count,'#04x').replace('0x','')
                product = plfam_dict[gid][plfam_id]['product']
                plfam_str = f'{plfam_id}\t{gid}\t{feature_count}\t{genome_count}\t{product}\t{aa_length_min}\t{aa_length_max}\t{aa_length_mean}\t{aa_length_std}\t{genomes}'
                plfam_line_list.append(plfam_str)
    
    for pgfam_id in pgfam_dict['unique_set']: 
        for gid in genome_ids: 
            if gid not in pgfam_dict:
                continue
            pgfam_data = {}
            pgfam_data['pgfam_id'] = pgfam_id
            pgfam_data['genome_id'] = gid
            if pgfam_id in pgfam_dict[gid]:
                aa_length_list = pgfam_dict[gid][pgfam_id]['aa_length_list'] 
                aa_length_max = max(aa_length_list)
                aa_length_min = min(aa_length_list)
                aa_length_mean = np.mean(aa_length_list)
                aa_length_std = np.std(aa_length_list)
                feature_count = pgfam_dict[gid][pgfam_id]['feature_count']
                genome_count = pgfam_dict[gid][pgfam_id]['genome_count']
                genomes = format(feature_count,'#04x').replace('0x','')
                product = pgfam_dict[gid][pgfam_id]['product']
                pgfam_str = f'{pgfam_id}\t{gid}\t{feature_count}\t{genome_count}\t{product}\t{aa_length_min}\t{aa_length_max}\t{aa_length_mean}\t{aa_length_std}\t{genomes}'
                pgfam_line_list.append(pgfam_str)
    
    output_json = {}
    output_json['plfam'] = '\n'.join(plfam_line_list) 
    output_json['pgfam'] = '\n'.join(pgfam_line_list) 
    #output_json['genome_ids'] = genome_ids
    output_json['genome_ids'] = list(set(genome_ids).intersection(present_genome_ids)) 
    tmp_data = genome_data.loc[genome_data['Genome ID'].isin(output_json['genome_ids'])]
    tmp_data.set_index('Genome ID',inplace=True)
    tmp_data = tmp_data.loc[output_json['genome_ids']]
    output_json['genome_names'] = tmp_data['Genome Name'].tolist()
    output_json['job_name'] = output_file

    output_json_file = os.path.join(output_dir,output_file+'_proteinfams_tables.json')
    with open(output_json_file,"w") as o:
        o.write(json.dumps(output_json))

    print("ProteinFamilies Complete")

def run_subsystems(genome_ids, query_dict, output_file, output_dir, genome_data, session):
    
    # json(facet,{"stat":{"type":"field","field":"superclass","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","class":{"type":"field","field":"class","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","gene_count":"unique(feature_id)","subclass":{"type":"field","field":"subclass","limit":-1,"facet":{"subsystem_count":"unique(subsystem_id)","gene_count":"unique(feature_id)"}}}}}}}):  

    # subsystems_df = getSubsystemsDataFrame(genome_ids,session) 
    if 'subsystems' not in query_dict:
        return ({ 'success': False })
    subsystems_df = query_dict['subsystems']
    
    # Superclass, class, and subclass can be different cases: convert all to lower case
    #subsystems_df['superclass'] = subsystems_df['superclass'].str.lower()
    #subsystems_df['class'] = subsystems_df['class'].str.lower()
    #subsystems_df['subclass'] = subsystems_df['subclass'].str.lower()

    # query 
    #base_query = "https://www.patricbrc.org/api/subsystem/?in(genome_id,("
    #end_query = "))&limit(200000000)&http_accept=text/tsv"
    #query = base_query + ",".join(genome_ids) + end_query
    #print("Subsystems Query:\n{0}".format(query))
    #subsystems_df = pd.read_csv(query,sep="\t")
    subsystems_file = os.path.join(output_dir,output_file+"_subsystems.tsv")
    subsystems_df.to_csv(subsystems_file, header=True, sep="\t")

    # TODO: remove, used for testing
    #subsystems_df = pd.read_csv(subsystems_file, sep="\t", index_col=0)

    # faceting for subsystems overview
    overview_dict = {}
    key_set = set()
    subsystem_present_genomes = set()
    for genome_id in subsystems_df['genome_id'].unique():
        subsystem_present_genomes.add(genome_id)
        genome_df = subsystems_df.loc[subsystems_df['genome_id'] == genome_id]
        overview_dict[genome_id] = {}
        overview_dict[genome_id]["superclass_counts"] = len(genome_df['subsystem_id'].unique())
        # TODO: check that this is correct for each level
        overview_dict[genome_id]["gene_counts"] = genome_df.shape[0] 
        for superclass in genome_df['superclass'].str.lower().unique():
            if type(superclass) != str:
                continue
            superclass_df = subsystems_df.loc[subsystems_df['superclass'].str.lower() == superclass]
            overview_dict[genome_id][superclass] = {}
            overview_dict[genome_id][superclass]["class_counts"] = len(superclass_df['subsystem_id'].unique())
            overview_dict[genome_id][superclass]["gene_counts"] = superclass_df.shape[0] 
            for clss in superclass_df['class'].str.lower().unique():
                if type(clss) != str:
                    continue
                class_df = superclass_df.loc[superclass_df['class'].str.lower() == clss]
                overview_dict[genome_id][superclass][clss] = {}
                overview_dict[genome_id][superclass][clss]['subclass_counts'] = len(class_df['subsystem_id'].unique())
                overview_dict[genome_id][superclass][clss]['gene_counts'] = class_df.shape[0] 
                for subclass in class_df['subclass'].str.lower().unique():
                    if type(subclass) != str:
                        continue
                    subclass_df = class_df.loc[class_df['subclass'].str.lower() == subclass]
                    overview_dict[genome_id][superclass][clss][subclass] = {}
                    overview_dict[genome_id][superclass][clss][subclass]['subsystem_name_counts'] = len(subclass_df['subsystem_id'].unique()) 
                    overview_dict[genome_id][superclass][clss][subclass]['gene_counts'] = subclass_df.shape[0] 
                    #for name in subsystems_df['subsystem_name'].unique(): 
                    #    key = build_key([superclass,clss,subclass,name]) #TODO: might not need 
                    #    key_set.add(key)

    #subsystems_overview_file = subsystems_file.replace('.tsv','_overview.json')
    #with open(subsystems_overview_file,'w') as o:
    #    json.dump(overview_dict,o)

    # faceting for subsystems table
    st_list = [] #subsystem table list
    # get stats on a per-genome basis
    # use unique subsystems_ids to create subsystem table for each genome
    for genome_id in subsystems_df['genome_id'].unique(): 
        print('---{0}'.format(genome_id))
        genome_df = subsystems_df.loc[subsystems_df['genome_id'] == genome_id]
        #genome_table = genome_df.drop(['feature_id','public','refseq_locus_tag','role_id','genome_id','taxon_id','gene','role_name','owner','product','patric_id','genome_name','id','_version_','date_inserted','date_modified'], axis=1)
        genome_table = genome_df.drop(return_columns_to_remove('subsystems_subsystems',genome_df.columns.tolist()), axis=1)
        genome_table = genome_table.drop_duplicates()
        
        genome_table['genome_count'] = [1]*genome_table.shape[0]
        genome_table['gene_count'] = [0]*genome_table.shape[0]
        genome_table['role_count'] = [0]*genome_table.shape[0]
    
        genome_table['genome_id'] = [genome_id]*genome_table.shape[0]

        # TODO: unknown if this needs to stay or go
        # Get unique subsystem ids, first row for information
        keep_rows = []
        sub_id_list = []
        for i in range(0,genome_table.shape[0]):
            if not genome_table.iloc[i]['subsystem_id'] in sub_id_list:
                sub_id_list.append(genome_table.iloc[i]['subsystem_id'])
                keep_rows.append(i)
        genome_table = genome_table.iloc[keep_rows]

        # add stats columns 
        for sub_id in genome_table['subsystem_id']:
            tmp_df = genome_df.loc[sub_id] 
            if isinstance(tmp_df, pd.DataFrame): # if only one record, returns as Series
                if 'gene' in tmp_df.columns.tolist():
                    genome_table.loc[sub_id,'gene_count'] = len(tmp_df['gene'].unique()) # unsure if this is correct)
                else: # unsure if this is correct
                    genome_table.loc[sub_id,'gene_count'] = 0 
                genome_table.loc[sub_id,'role_count'] = len(tmp_df['role_id'].unique())
            else:
                genome_table.loc[sub_id,'gene_count'] = 1
                genome_table.loc[sub_id,'role_count'] = 1
        
        st_list.append(genome_table)    

    subsystems_table = pd.concat(st_list)
    #subsystems_table_output_file = subsystems_file.replace('.tsv','_summary.tsv')
    #subsystems_table.to_csv(subsystems_table_output_file,sep="\t",index=False)

    gene_df = query_dict['feature']
    # change column names
    
    #gene_df.drop(['_version_'],inplace=True)
    # Add subsystems columns to genes table
    #remove_columns = ['date_inserted','date_modified','genome_name','gene','owner','patric_id','public','product','refseq_locus_tag','taxon_id','_version_']
    gene_df = pd.merge(gene_df,subsystems_df.drop(return_columns_to_remove('subsystems_genes',subsystems_df.columns.tolist()),axis=1),on=['genome_id','feature_id'],how='inner')

    output_json_file = subsystems_file.replace('.tsv','_tables.json')
    
    output_json = {}
    output_json['genome_ids'] = list(subsystem_present_genomes) 
    output_json['genome_names'] = genome_data.set_index('Genome ID').loc[list(subsystem_present_genomes)]['Genome Name'].tolist() # returns a list of genome names in the same order as the genome ids 
    output_json['overview'] = overview_dict 
    output_json['job_name'] = output_file
    output_json['subsystems'] = subsystems_table.to_csv(index=False,sep='\t')
    output_json['genes'] = gene_df.to_csv(index=False,sep='\t')
    with open(output_json_file,'w') as o:
        o.write(json.dumps(output_json))

    print('Subsystems complete')
    return ({ 'success': True, 'genomes': list(subsystem_present_genomes) })

def run_subsystems_v2(genome_ids, query_dict, output_file, output_dir, genome_data, session):

    subsystems_file = os.path.join(output_dir,output_file+'_subsystems.tsv')
    subsystem_line_list = []
    subsystem_header = 'superclass\tclass\tsubclass\tsubsystem_name\tgene_count\trole_count'
    subsystem_line_list.append(subsystem_header)
    subsystem_dict = {}
    overview_counts_dict = {}

    subsystem_query_data = []
    required_fields = ['superclass','class','subclass','subsystem_name','feature_id','gene','product','role_id','role_name']
    subsystem_data_found = False
    subsystem_genomes_found = set()
    subsystem_table_header = None
    print_one = True
    for gids in chunker(genome_ids, 20):
        base = "https://www.patricbrc.org/api/subsystem/?http_download=true"
        query = f"in(genome_id,({','.join(gids)}))&limit(2500000)&sort(+id)"
        headers = {"accept":"application/json", "content-type":"application/rqlquery+x-www-form-urlencoded","Authorization": session.headers['Authorization']}

        #dict_keys(['active', 'class', 'date_inserted', 'date_modified', 'feature_id', 'gene', 'genome_id', 'genome_name', 'id', 'owner', 'patric_id', 'product', 'public', 'refseq_locus_tag', 'role_id', 'role_name', 'subclass', 'subsystem_id', 'subsystem_name', 'superclass', 'taxon_id', '_version_'])

        print('Query = {0}\nHeaders = {1}'.format(base+'&'+query,headers))
        result_header = True        
        current_header = None
        all_data = json.loads(getQueryDataText(base,query,headers))
        for line in all_data:
            subsystem_data_found = True
            if result_header:
                result_header = False
                print(line)
                current_header = line.keys()
                if subsystem_table_header is None or len(current_header) > len(subsystem_table_header):
                    subsystem_table_header = current_header
                continue
            if print_one:
                print_one = False
                print(line)
            subsystem_fields = line
            for field in required_fields:
                if field not in subsystem_fields:
                    subsystem_fields[field] = ''
            subsystem_query_data.append(subsystem_fields)
            try:
                active = subsystem_fields['active'] 
                clss = subsystem_fields['class'] 
                feature_id = subsystem_fields['feature_id'] 
                gene = subsystem_fields['gene'] 
                genome_id = subsystem_fields['genome_id']
                product = subsystem_fields['product'] 
                role_id = subsystem_fields['role_id'] 
                role_name = subsystem_fields['role_name'] 
                subclass = subsystem_fields['subclass'] 
                subsystem_id = subsystem_fields['subsystem_id'] 
                subsystem_name = subsystem_fields['subsystem_name'] 
                superclass = subsystem_fields['superclass']
                # TODO: underlying issue of metadata, capitalizition of superclasses is not consistent
                superclass = superclass.upper()
            except Exception as e:
                sys.stderr.write(f'Error with the following line:\n{e}\n{line}\n')
                continue
            subsystem_genomes_found.add(genome_id)
            # TODO: chat with Maulik about the correct gene/role count
            # - Do I need to keep track of unique features and their counts???
            if superclass not in subsystem_dict:
                subsystem_dict[superclass] = {} 
                overview_counts_dict[superclass] = {}
            if clss not in subsystem_dict[superclass]:
                subsystem_dict[superclass][clss] = {}
                overview_counts_dict[superclass][clss] = {}
            if subclass not in subsystem_dict[superclass][clss]:
                subsystem_dict[superclass][clss][subclass] = {}
                overview_counts_dict[superclass][clss][subclass] = {}
                overview_counts_dict[superclass][clss][subclass]['subsystem_names'] = set()
                overview_counts_dict[superclass][clss][subclass]['gene_set'] = set()
            if subsystem_name not in subsystem_dict[superclass][clss][subclass]:
                subsystem_dict[superclass][clss][subclass][subsystem_name] = {}
                subsystem_dict[superclass][clss][subclass][subsystem_name]['gene_set'] = set()
                subsystem_dict[superclass][clss][subclass][subsystem_name]['role_set'] = set()
            overview_counts_dict[superclass][clss][subclass]['subsystem_names'].add(subsystem_name)
            if gene != '' or gene is not None:
                subsystem_dict[superclass][clss][subclass][subsystem_name]['gene_set'].add(gene)
                overview_counts_dict[superclass][clss][subclass]['gene_set'].add(gene)
            if role_id != '' or gene is not None: 
                subsystem_dict[superclass][clss][subclass][subsystem_name]['role_set'].add(role_id)

    # gets counts for overview dict, any other adjustments
    # create subsystems table
    subsystems_table_list = []
    overview_dict = {} 
    for superclass in overview_counts_dict:
        overview_dict[superclass] = {}
        overview_dict[superclass]['subsystem_name_counts'] = 0
        overview_dict[superclass]['gene_counts'] = 0
        for clss in overview_counts_dict[superclass]:
            overview_dict[superclass][clss] = {}
            overview_dict[superclass][clss]['subsystem_name_counts'] = 0
            overview_dict[superclass][clss]['gene_counts'] = 0
            for subclass in overview_counts_dict[superclass][clss]:
                overview_dict[superclass][clss][subclass] = {}
                overview_dict[superclass][clss][subclass]['subsystem_name_counts'] = len(overview_counts_dict[superclass][clss][subclass]['subsystem_names'])
                overview_dict[superclass][clss][subclass]['gene_counts'] = len(overview_counts_dict[superclass][clss][subclass]['gene_set'])
                overview_dict[superclass][clss]['subsystem_name_counts'] += len(overview_counts_dict[superclass][clss][subclass]['subsystem_names'])
                overview_dict[superclass][clss]['gene_counts'] += len(overview_counts_dict[superclass][clss][subclass]['gene_set'])
                overview_dict[superclass]['gene_counts'] += len(overview_counts_dict[superclass][clss][subclass]['gene_set'])
                overview_dict[superclass]['subsystem_name_counts'] += len(overview_counts_dict[superclass][clss][subclass]['subsystem_names'])
                for subsystem_name in subsystem_dict[superclass][clss][subclass]:
                    new_entry = {
                        'superclass': superclass,
                        'class': clss,
                        'subclass': subclass,
                        'subsystem_name': subsystem_name,
                        'role_counts': len(subsystem_dict[superclass][clss][subclass][subsystem_name]['role_set']),
                        'gene_counts': len(subsystem_dict[superclass][clss][subclass][subsystem_name]['gene_set'])
                    }
                    subsystems_table_list.append(new_entry)

    if not subsystem_data_found:
        return ({ 'success': False }) 

    parsed_query_data = []
    for line in subsystem_query_data:
        new_line = ''
        for field in subsystem_table_header:
            if new_line != '':
                new_line += '\t'
            if field not in line:    
                new_line += ' '
            else:
                value = line[field]
                if not isinstance(value,str):
                    value = str(value)
                new_line += value 
        parsed_query_data.append(new_line.split('\t'))

    subsystem_df = pd.DataFrame(parsed_query_data,columns=subsystem_table_header)
    subsystems_table = pd.DataFrame(subsystems_table_list)

    gene_df = query_dict['feature']
    gene_df = pd.merge(gene_df,subsystem_df.drop(return_columns_to_remove('subsystems_genes',subsystem_df.columns.tolist()),axis=1),on=['genome_id','feature_id'],how='inner')
    
    output_json_file = subsystems_file.replace('.tsv','_tables.json')
    
    output_json = {}
    output_json['genome_ids'] = list(subsystem_genomes_found)
    output_json['genome_names'] = genome_data.set_index('Genome ID').loc[list(subsystem_genomes_found)]['Genome Name'].tolist() # returns a list of genome names in the same order as the genome ids
    output_json['overview'] = overview_dict
    output_json['job_name'] = output_file
    output_json['subsystems'] = subsystems_table.to_csv(index=False,sep='\t')
    output_json['genes'] = gene_df.to_csv(index=False,sep='\t')
    with open(output_json_file,'w') as o:
        o.write(json.dumps(output_json))

    print('Subsystems complete')
    return ({ 'success': True, 'genomes': list(subsystem_genomes_found) })

def run_pathways_v2(genome_ids, query_dict, output_file, output_dir, genome_data, session):
    
    pathways_file = os.path.join(output_dir,output_file+'_pathways.tsv')
    #pathway_df = query_dict['pathway']
    #pathway_df.to_csv(pathways_file,sep='\t',index=False)
    pathway_line_list = []
    ec_line_list = []
    pathway_header = 'annotation\tpathway_id\tpathway_name\tpathway_class\tgenome_count\tec_count\tgene_count\tgenome_ec\tec_conservation\tgene_conservation'
    ec_header = 'annotation\tpathway_id\tpathway_name\tpathway_class\tec_description\tec_number\tgenome_count\tec_count\tgene_count\tgenome_ec'
    pathway_line_list.append(pathway_header)
    ec_line_list.append(ec_header)
    pathway_dict = {}
    ec_dict = {}
    unique_pathways = set()
    unique_ecs = set()
    unique_features = set()
    unique_pathway_features = {} 
    unique_pathway_ecs = {}
    
    pathway_query_data = []
    required_fields = ['annotation','ec_description','ec_number','feature_id','genome_id','pathway_class','pathway_id','pathway_name','patric_id','product']
    pathway_data_found = False
    pathway_genomes_found = set()
    pathway_table_header = None
    for gids in chunker(genome_ids, 20):
        base = "https://www.patricbrc.org/api/pathway/?http_download=true"
        query = f"in(genome_id,({','.join(gids)}))&limit(2500000)&sort(+id)&eq(annotation,PATRIC)"
        headers = {"accept":"application/json", "content-type":"application/rqlquery+x-www-form-urlencoded", "Authorization": session.headers['Authorization']}
        
        print('Query = {0}\nHeaders = {1}'.format(base+'&'+query,headers))
        #accession       alt_locus_tag   annotation      date_inserted   date_modified   ec_description  ec_number       feature_id      genome_ec       genome_id       genome_name     id      owner   pathway_class   pathway_ec      pathway_id   pathway_name     patric_id       product public  refseq_locus_tag        sequence_id     taxon_id        _version_

        result_header = True
        current_header = None
        all_data = json.loads(getQueryDataText(base,query,headers))
        for line in all_data:
            pathway_data_found = True
            if result_header:
                result_header = False
                print(line)
                current_header = line.keys()
                if pathway_table_header is None or len(current_header) > len(pathway_table_header):
                    pathway_table_header = current_header 
                #pathway_query_data.append(line)
                continue
            #line = line.strip().replace('\"','').split('\t')
            pathway_fields = line
            #for idx,f in enumerate(line):
            #    pathway_fields[pathway_table_header[idx]] = f
            for field in required_fields:
                if field not in pathway_fields:
                    pathway_fields[field] = ''
            pathway_query_data.append(pathway_fields)
            try:
                annotation = pathway_fields['annotation'] 
                ec_description = pathway_fields['ec_description'] 
                ec_number = pathway_fields['ec_number'] 
                feature_id = pathway_fields['feature_id'] 
                genome_id = pathway_fields['genome_id'] 
                #genome_name = line[10].replace('\"','')
                pathway_class =  pathway_fields['pathway_class'] 
                pathway_id = pathway_fields['pathway_id'] 
                pathway_name = pathway_fields['pathway_name'] 
                patric_id = pathway_fields['patric_id'] 
                product = pathway_fields['product'] 
            except Exception as e:
                sys.stderr.write(f'Error with the following line:\n{e}\n{line}\n')
                continue
            pathway_genomes_found.add(genome_id)
            unique_pathways.add(pathway_id)
            unique_ecs.add(ec_number)
            #unique_features.add(patric_id)
            '''
            if pathway_id not in unique_pathway_features:
                unique_pathway_features[pathway_id] = {} 
            if 'gene' in pathway_fields:
                pathway_gene = pathway_fields['gene']
                if pathway_gene not in unique_pathway_features[pathway_id]: 
                    unique_pathway_features[pathway_id][pathway_gene] = set()
                unique_pathway_features[pathway_id][pathway_gene].add(genome_id)
            '''
            if pathway_id not in unique_pathway_ecs:
                unique_pathway_ecs[pathway_id] = {}
            if ec_number not in unique_pathway_ecs[pathway_id]:
                unique_pathway_ecs[pathway_id][ec_number] = set()
            unique_pathway_ecs[pathway_id][ec_number].add(genome_id)

            # pathway data
            if pathway_id not in pathway_dict:
                pathway_dict[pathway_id] = {} 
                pathway_dict[pathway_id]['annotation'] = annotation 
                pathway_dict[pathway_id]['pathway_id'] = pathway_id
                pathway_dict[pathway_id]['pathway_name'] = pathway_name
                pathway_dict[pathway_id]['pathway_class'] = pathway_class
                pathway_dict[pathway_id]['genome_count'] = set()
                pathway_dict[pathway_id]['ec_count'] = set()
                pathway_dict[pathway_id]['gene_count'] = set()
                pathway_dict[pathway_id]['genome_ec'] = set() 
            pathway_dict[pathway_id]['genome_count'].add(genome_id)
            pathway_dict[pathway_id]['ec_count'].add(ec_number)
            pathway_dict[pathway_id]['gene_count'].add(feature_id)
            pathway_dict[pathway_id]['genome_ec'].add(genome_id+'_'+ec_number)
            # ec data
            #ec_header = 'annotation\tpathway_id\tpathway_name\tpathway_class\tproduct\tec_number\tgenome_count\tec_count\tgene_count\tgenome_ec'
            if pathway_id not in ec_dict:
                ec_dict[pathway_id] = {}
            if ec_number not in ec_dict[pathway_id]:
                ec_dict[pathway_id][ec_number] = {}
                ec_dict[pathway_id][ec_number]['annotation'] = annotation
                ec_dict[pathway_id][ec_number]['pathway_id'] = pathway_id
                ec_dict[pathway_id][ec_number]['pathway_name'] = pathway_name
                ec_dict[pathway_id][ec_number]['pathway_class'] = pathway_class
                ec_dict[pathway_id][ec_number]['ec_description'] = ec_description 
                ec_dict[pathway_id][ec_number]['ec_number'] = ec_number
                ec_dict[pathway_id][ec_number]['genome_count'] = set()
                ec_dict[pathway_id][ec_number]['ec_count'] = set()
                ec_dict[pathway_id][ec_number]['gene_count'] = set()
                ec_dict[pathway_id][ec_number]['genome_ec'] = set()
            ec_dict[pathway_id][ec_number]['genome_count'].add(genome_id)
            ec_dict[pathway_id][ec_number]['ec_count'].add(ec_number)
            ec_dict[pathway_id][ec_number]['gene_count'].add(feature_id)
            ec_dict[pathway_id][ec_number]['genome_ec'].add(genome_id+'_'+ec_number)

    if not pathway_data_found:
        return ({ 'success': False }) 

    parsed_query_data = []
    for line in pathway_query_data:
        new_line = ''
        for field in pathway_table_header:
            if new_line != '':
                new_line += '\t'
            if field not in line:    
                new_line += ' '
            else:
                value = line[field]
                if not isinstance(value,str):
                    value = str(value)
                new_line += value 
        parsed_query_data.append(new_line.split('\t'))

    pathway_df = pd.DataFrame(parsed_query_data,columns=pathway_table_header)
    gene_df = query_dict['feature']

    genes_output = pd.merge(gene_df.drop(return_columns_to_remove('pathways_genes',gene_df.columns.tolist()), axis=1),pathway_df,on=['genome_id','patric_id'],how='inner')


    if 'gene_x' in genes_output.columns:
        genes_output['gene'] = genes_output['gene_x']
        genes_output.drop(['gene_x','gene_y'],inplace=True,axis=1)

    for idx in range(0,genes_output.shape[0]):
        pathway_id = genes_output.iloc[idx].pathway_id
        if pathway_id not in unique_pathway_features:
            unique_pathway_features[pathway_id] = {}    
        gene = genes_output.iloc[idx]['gene']
        if gene is None or gene is np.nan:
            continue
        genome_id = genes_output.iloc[idx].genome_id
        if gene not in unique_pathway_features[pathway_id]:
            unique_pathway_features[pathway_id][gene] = set()
        unique_pathway_features[pathway_id][gene].add(genome_id)
        unique_features.add(gene)

    # get gene data frame 
    # get conservation stats and add lines
    for pathway_id in pathway_dict:
        pathway_dict[pathway_id]['genome_count'] = len(pathway_dict[pathway_id]['genome_count'])
        pathway_dict[pathway_id]['ec_count'] = len(pathway_dict[pathway_id]['ec_count'])
        pathway_dict[pathway_id]['gene_count'] = len(pathway_dict[pathway_id]['gene_count'])
        pathway_dict[pathway_id]['genome_ec'] = len(pathway_dict[pathway_id]['genome_ec'])
        #pathway_dict[pathway_id]['ec_conservation'] = float(len(unique_pathway_ecs[pathway_id]))/float(len(unique_ecs))*100.0
        #pathway_dict[pathway_id]['gene_conservation'] = float(len(unique_pathway_features[pathway_id]))/float(len(unique_features))*100.0
        annotation = pathway_dict[pathway_id]['annotation']
        pathway_id = pathway_dict[pathway_id]['pathway_id']
        pathway_name = pathway_dict[pathway_id]['pathway_name']
        pathway_class = pathway_dict[pathway_id]['pathway_class']
        genome_count = pathway_dict[pathway_id]['genome_count']
        ec_count = pathway_dict[pathway_id]['ec_count']
        gene_count = pathway_dict[pathway_id]['gene_count']
        genome_ec = pathway_dict[pathway_id]['genome_ec']
        # calculate ec_conservation score
        ec_numerator = 0
        ec_denominator = 0
        for ec_number in unique_pathway_ecs[pathway_id]:
            ec_numerator += len(unique_pathway_ecs[pathway_id][ec_number])
            ec_denominator += len(pathway_genomes_found)
        #ec_numerator = float(ec_numerator) * float(len(pathway_genomes_found))
        #ec_denominator = float(ec_denominator) * float(len(pathway_genomes_found))
        if ec_denominator == 0:
            ec_conservation = 0
        else:
            ec_conservation = float(ec_numerator) / float(ec_denominator) * 100.0
        # calculate gene_conservation
        gene_numerator = 0
        gene_denominator = 0
        for gene in unique_pathway_features[pathway_id]:
            gene_numerator += len(unique_pathway_features[pathway_id][gene])
            gene_denominator += len(pathway_genomes_found)
        if gene_denominator == 0:
            gene_conservation = 0
        else:
            gene_conservation = float(gene_numerator) / float(gene_denominator) * 100.0
        pathway_line = f'{annotation}\t{pathway_id}\t{pathway_name}\t{pathway_class}\t{genome_count}\t{ec_count}\t{gene_count}\t{genome_ec}\t{ec_conservation}\t{gene_conservation}'
        pathway_line_list.append(pathway_line)
        # now EC data
        for ec_number in ec_dict[pathway_id]:
            ec_dict[pathway_id][ec_number]['genome_count'] = len(ec_dict[pathway_id][ec_number]['genome_count'])
            ec_dict[pathway_id][ec_number]['ec_count'] = len(ec_dict[pathway_id][ec_number]['ec_count'])
            ec_dict[pathway_id][ec_number]['gene_count'] = len(ec_dict[pathway_id][ec_number]['gene_count'])
            ec_dict[pathway_id][ec_number]['genome_ec'] = len(ec_dict[pathway_id][ec_number]['genome_ec'])
            annotation = ec_dict[pathway_id][ec_number]['annotation']
            pathway_id = ec_dict[pathway_id][ec_number]['pathway_id']
            pathway_name = ec_dict[pathway_id][ec_number]['pathway_name']
            pathway_class = ec_dict[pathway_id][ec_number]['pathway_class']
            ec_description = ec_dict[pathway_id][ec_number]['ec_description']
            ec_number = ec_dict[pathway_id][ec_number]['ec_number']
            genome_count = ec_dict[pathway_id][ec_number]['genome_count']
            ec_count = ec_dict[pathway_id][ec_number]['ec_count']
            gene_count = ec_dict[pathway_id][ec_number]['gene_count']
            genome_ec = ec_dict[pathway_id][ec_number]['genome_ec']
            ec_line = f'{annotation}\t{pathway_id}\t{pathway_name}\t{pathway_class}\t{ec_description}\t{ec_number}\t{genome_count}\t{ec_count}\t{gene_count}\t{genome_ec}'
            ec_line_list.append(ec_line)

    pathway_output = '\n'.join(pathway_line_list)
    ec_output = '\n'.join(ec_line_list)

    output_json = {}
    output_json['pathway'] = pathway_output
    output_json['ecnumber'] = ec_output
    output_json['genes'] = genes_output.to_csv(index=False,sep='\t')
    output_json['genome_ids'] = list(pathway_genomes_found) 
    output_json['job_name'] = output_file
    
    pathway_df.to_csv(pathways_file,sep='\t',index=False)

    output_json_file = pathways_file.replace('.tsv','_tables.json')
    with open(output_json_file,"w") as o:
        o.write(json.dumps(output_json))

    print("Pathways Complete")
    pathway_success_json = {
        'genomes': list(pathway_genomes_found),
        'success': True
    }
    return pathway_success_json

def run_pathways(genome_ids, query_dict, output_file,output_dir, genome_data, session):
    
    pathways_file = os.path.join(output_dir,output_file+'_pathways.tsv')
    # TODO: create alt_locus_tag query
    # pathway_df = getPathwayDataFrame(genome_ids,session, limit=2500000)
    pathway_df = query_dict['pathway']
    pathway_df.to_csv(pathways_file,sep='\t',index=False)

    #TODO: remove, reading in file for testing
    #pathway_df = pd.read_csv(pathways_file,sep="\t")
    pathways_list = []
    ecnum_list = []
    genes_info_dict = {} # key is patric_id
    ec_con_table = {}
    # TODO:
    # - make sure to check genome_id type issue 
    for genome_id in genome_ids: 
        print('---Faceting GenomeId: {0}---'.format(genome_id))
        genome_df = pathway_df.loc[pathway_df['genome_id'] == genome_id]

        '''
        pathway_table = genome_df.drop(['pathway_ec','genome_name','accession','genome_ec',
                                        'product','gene','public','patric_id','sequence_id','ec_number',
                                        'feature_id','taxon_id','ec_description','refseq_locus_tag','owner',
                                        'id','_version_','date_inserted','date_modified'], axis=1)
        # TODO: add index column
        ec_table = genome_df.drop(['pathway_ec','genome_name','accession','genome_ec','product','feature_id',
                                    'gene','public','patric_id','sequence_id','taxon_id','refseq_locus_tag',
                                    'owner','id','_version_','date_inserted','date_modified'], axis=1)
        '''
        pathway_table = genome_df[['genome_id','annotation','pathway_class','pathway_name','pathway_id']]
        ec_table = genome_df[['genome_id','annotation','pathway_class','pathway_name','pathway_id','ec_number','ec_description','ec_index']]

        pathway_table = pathway_table.drop_duplicates()
        ec_table = ec_table.drop_duplicates()

        # add pathway stats columns 
        pathway_table['genome_count'] = [1]*pathway_table.shape[0]
        pathway_table['gene_count'] = [0]*pathway_table.shape[0]
        pathway_table['ec_count'] = [0]*pathway_table.shape[0]
        pathway_table['genome_ec'] = [0]*pathway_table.shape[0]
        for pathway_id in pathway_table['pathway_id']: # should be unique
            tmp_df = genome_df.loc[pathway_id]
            if isinstance(tmp_df, pd.DataFrame): # if only one entry, returns a Series 
                pathway_table.loc[pathway_id,'gene_count'] = len(tmp_df['feature_id'].unique())
                pathway_table.loc[pathway_id,'ec_count'] = len(tmp_df['ec_number'].unique())
                pathway_table.loc[pathway_id,'genome_ec'] = len(tmp_df['genome_ec'].unique())
                # for genes info: take first record
                p_id = tmp_df.iloc[0]['patric_id']
            else:
                pathway_table.loc[pathway_id,'gene_count'] = 1 
                pathway_table.loc[pathway_id,'ec_count'] = 1 
                pathway_table.loc[pathway_id,'genome_ec'] = 1 
                # for genes info: take first record
                p_id = tmp_df['patric_id']
            if p_id not in genes_info_dict:
                genes_info_dict[p_id] = tmp_df.iloc[0] 

        # TODO: maybe modify once the multiple-rows thing has been decided
        # get first ec_number entry data for ec_number duplicates
        # pathway_id is not duplicated for each ec_number
        keep_rows = []
        ec_list = []
        for i in range(0,ec_table.shape[0]):
            if not ec_table.iloc[i]['ec_number'] in ec_list:
                ec_list.append(ec_table.iloc[i]['ec_number'])
                keep_rows.append(i)
        ec_table = ec_table.iloc[keep_rows]

        # add ec_number stats columns
        ec_table['genome_count'] = [1]*ec_table.shape[0]
        ec_table['gene_count'] = [0]*ec_table.shape[0]
        ec_table['ec_count'] = [0]*ec_table.shape[0]
        ec_table['genome_ec'] = [0]*ec_table.shape[0]
        genome_df.set_index('ec_index',inplace=True)
        ec_table.set_index('ec_index',inplace=True)
        for ec_number in ec_table['ec_number']:
            tmp_df = genome_df.loc[ec_number]
            if isinstance(tmp_df, pd.DataFrame): # if only one entry, returns a Series
                ec_table.loc[ec_number,'gene_count'] = len(tmp_df['feature_id'].unique()) 
                ec_table.loc[ec_number,'ec_count'] = len(tmp_df['ec_number'].unique()) 
                ec_table.loc[ec_number,'genome_ec'] = len(tmp_df['ec_number'].unique())
            else:
                ec_table.loc[ec_number,'gene_count'] = 1 
                ec_table.loc[ec_number,'ec_count'] = 1 
                ec_table.loc[ec_number,'genome_ec'] = 1 
        # append to lists
        pathways_list.append(pathway_table)
        ecnum_list.append(ec_table)

    # concatenate tables and do final calculations: genome_count, genome_ec
    # counting is done per-genome, multi-genome calculation adjustments are done on the front end
    pathway_output = pd.concat(pathways_list)
    ec_output = pd.concat(ecnum_list)

    # Get gene data
    feature_list = []
    # gene_df = getFeatureDataFrame(genome_ids,session, limit=2500000)
    gene_df = query_dict['feature']
    
    genes_output = pd.merge(gene_df.drop(return_columns_to_remove('pathways_genes',gene_df.columns.tolist()), axis=1),pathway_df,on=['genome_id','patric_id'],how='inner')

    if False:
        # Parse gene data
        # TODO:
        # - make sure to check genome_id type issue 
        genes_list = []
        for genome_id in genome_ids:
            print('---Faceting GenomeId Genes Table: {0}---'.format(genome_id))
            genome_df = gene_df.loc[gene_df['genome_id'] == genome_id]
            
            genes_table = genome_df.drop_duplicates()
            genes_table.gene = genes_table.gene.fillna('')
        
            # get first gene entry data for gene duplicates
            # pathway_id is not duplicated for each gene
            keep_rows = []
            g_list = []
            for i in range(0,genes_table.shape[0]):
                if not genes_table.iloc[i]['gene'] in g_list:
                    g_list.append(genes_table.iloc[i]['gene'])
                    keep_rows.append(i)
                if genes_table.iloc[i]['gene'] == '':
                    keep_rows.append(i)
            genes_table = genes_table.iloc[keep_rows]

            # fill in ec_description, ec_number, index, pathwayid, pathway name 
            genes_table['ec_description'] = ['']*genes_table.shape[0]
            genes_table['ec_number'] = ['']*genes_table.shape[0]
            genes_table['index'] = ['']*genes_table.shape[0]
            genes_table['pathway_id'] = ['']*genes_table.shape[0]
            genes_table['pathway_name'] = ['']*genes_table.shape[0]
            # genes table stats
            genes_table['genome_count'] = [1]*genes_table.shape[0]
            genes_table['gene_count'] = [0]*genes_table.shape[0]
            genes_table['ec_count'] = [0]*genes_table.shape[0]
            genes_table['genome_ec'] = [0]*genes_table.shape[0]
            #genes_table['alt_locus_tag'] = ['']*genes_table.shape[0]
            # use patric_id to account for '' genes
            for p_id in genes_table['patric_id']:
                tmp_df = pathway_df.loc[pathway_df['genome_id'] == genome_id]
                tmp_df = tmp_df[tmp_df['patric_id'] == p_id]
                genes_table.loc[genes_table['patric_id'] == p_id,'gene_count'] = len(tmp_df['feature_id'].unique())
                genes_table.loc[genes_table['patric_id'] == p_id,'ec_count'] = len(tmp_df['ec_number'].unique())
                genes_table.loc[genes_table['patric_id'] == p_id,'genome_ec'] = len(tmp_df['ec_number'].unique())
                if p_id in genes_info_dict:
                    genes_table.loc[genes_table['patric_id'] == p_id,'ec_description'] = genes_info_dict[p_id]['ec_description']
                    genes_table.loc[genes_table['patric_id'] == p_id,'ec_number'] = genes_info_dict[p_id]['ec_number']
                    genes_table.loc[genes_table['patric_id'] == p_id,'index'] = genes_info_dict[p_id]['id']
                    genes_table.loc[genes_table['patric_id'] == p_id,'pathway_id'] = genes_info_dict[p_id]['pathway_id']
                    genes_table.loc[genes_table['patric_id'] == p_id,'pathway_name'] = genes_info_dict[p_id]['pathway_name']
     
            genes_list.append(genes_table)

    #genes_output = pd.concat(genes_list)

    output_json = {}
    output_json['pathway'] = pathway_output.to_csv(index=False,sep='\t')
    output_json['ecnumber'] = ec_output.to_csv(index=False,sep='\t')
    output_json['genes'] = genes_output.to_csv(index=False,sep='\t')
    output_json['genome_ids'] = genome_ids
    output_json['job_name'] = output_file
    

    output_json_file = pathways_file.replace('.tsv','_tables.json')
    with open(output_json_file,"w") as o:
        o.write(json.dumps(output_json))

    print("Pathways Complete")

def generate_report(genome_ids, pathway_obj, subsystems_obj, proteinfams_obj, output_dir):
    report_text_list = []
    if pathway_obj['success']:
        report_text_list.append(f"Pathways succeded: {len(pathway_obj['genomes'])} out of {len(genome_ids)} genomes had pathway data") 
        if len(pathway_obj['genomes']) != len(genome_ids):
            missing_pathway_genomes = list(set(genome_ids).difference(set(pathway_obj['genomes'])))
            report_text_list.append(f"Genomes Missing from Pathways: {','.join(missing_pathway_genomes)}")
    else:
        report_text_list.append('Pathways Failed: see stdout and stderr')
    if subsystems_obj['success']:
        report_text_list.append(f"Subsystems succeeded: {len(subsystems_obj['genomes'])} out of {len(genome_ids)} genomes had subsystems data")
        if len(subsystems_obj['genomes']) != len(genome_ids):
            missing_subsystems_genomes = list(set(genome_ids).difference(set(subsystems_obj['genomes'])))
            report_text_list.append(f"Genomes Missing from Subsystems: {','.join(missing_subsystems_genomes)}")
    else:
        report_text_list.append('Subsystems Failed: see stdout and stderr')
    if proteinfams_obj['success']:
        report_text_list.append(f"ProteinFamilies succeeded: {len(proteinfams_obj['genomes'])} out of {len(genome_ids)} genomes had proteinfamilies data")
        if len(proteinfams_obj['genomes']) != len(genome_ids):
            missing_proteinfams_genomes = list(set(genome_ids).difference(set(proteinfams_obj['genomes'])))
            report_text_list.append(f"Genomes Missing from ProteinFamilies: {','.join(missing_proteinfams_genomes)}")
    else:
        report_text_list.append('ProteinFamilies Failed: see stdout and sterr')
    report_text = '\n'.join(report_text_list)
    report_file = os.path.join(output_dir,'report.txt')
    with open(report_file,'w') as o:
        o.write(report_text)
    

# Store pathways, subsystems, and features queries in a dictionary
def run_all_queries(genome_ids, session):
    query_dict = {}
    ### Run pathways query
    if False:
        print('pathways query')
        pathway_df = getPathwayDataFrame(genome_ids,session, limit=2500000)
        if not pathway_df is None:
            pathway_df['pathway_index'] = pathway_df['pathway_id']
            pathway_df['ec_index'] = pathway_df['ec_number']
            pathway_df.set_index('pathway_index', inplace=True)
            query_dict['pathway'] = pathway_df
        else:
            sys.stderr.write('Pathways dataframe is None\n')
    ### Run subsystems query
    if False:
        print('subsystems query')
        subsystems_df = getSubsystemsDataFrame(genome_ids,session) 
        if not subsystems_df is None:
            subsystems_df['subsystem_index'] = subsystems_df['subsystem_id']
            subsystems_df.set_index('subsystem_index', inplace=True)
            query_dict['subsystems'] = subsystems_df
        else:
            sys.stderr.write('Subsystems dataframe is None\n')
    ### Run features query
    if False:
        print('features query')
        feature_df = getFeatureDataFrame(genome_ids,session, limit=2500000)
        if not feature_df is None:
            column_map = {
                'Genome': 'genome_name',
                'Genome ID': 'genome_id',
                'Accession': 'accession',
                'BRC ID': 'patric_id',
                'RefSeq Locus Tag': 'refseq_locus_tag',
                'Alt Locus Tag': 'alt_locus_tag',
                'Feature ID': 'feature_id',
                'Annotation': 'annotation',
                'Feature Type': 'feature_type',
                'Start': 'start',
                'End': 'end',
                'Length': 'length',
                'Strand': 'strand',
                'FIGfam ID': 'figfam_id',
                'PATRIC genus-specific families (PLfams)': 'plfam_id',
                'PATRIC cross-genus families (PGfams)': 'pgfam_id',
                'Protein ID': 'protein_id',
                'AA Length': 'aa_length',
                'Gene Symbol': 'gene',
                'Product': 'product',
                'GO': 'go'
            }
            if 'Genome ID' in feature_df.columns:
                feature_df.rename(columns=column_map, inplace=True)
            feature_df['plfam_index'] = feature_df['plfam_id']
            feature_df['pgfam_index'] = feature_df['pgfam_id']
            #feature_df.set_index('plfam_index', inplace=True)
            query_dict['feature'] = feature_df
        else:
            sys.stderr.write('Features dataframe is None\n')
    return query_dict

def get_genome_group_ids(group_list,session):
    genome_group_ids = []
    for genome_group in group_list:
        genome_id_list = getGenomeIdsByGenomeGroup(genome_group,session,genomeGroupPath=True)
        genome_group_ids = genome_group_ids + genome_id_list
    return genome_group_ids

def run_compare_systems(job_data, output_dir):

    ###Setup session
    s = requests.Session()
    authenticateByEnv(s)

    output_dir = os.path.abspath(output_dir)
    if not os.path.exists(output_dir):
        subprocess.call(["mkdir", "-p", output_dir])
    output_file = job_data["output_file"]


    print("Run ComparativeSystems:\njob_data = {0}".format(job_data)) 
    print("output_dir = {0}".format(output_dir)) 
    
    genome_ids = job_data["genome_ids"]
    if len(job_data["genome_groups"]) > 0:
        genome_group_ids = get_genome_group_ids(job_data["genome_groups"],s)
        if len(genome_group_ids) == 0:
            sys.stderr.write('FAILED to get genome ids for genome groups: exiting')
            sys.exit(-1)
        # make ids unique 
        genome_ids = genome_ids + genome_group_ids
        genome_ids = list(set(genome_ids))

    # optionally add more genome info to output 
    genome_data = getDataForGenomes(genome_ids,s) 
    
    query_dict = run_all_queries(genome_ids, s)

    # TODO: add chunking
    # TODO: add recipe
    # TODO: add multithreading
    #pathway_success = run_pathways_v2(genome_ids, query_dict, output_file, output_dir, genome_data, s)
    #subsystems_success = run_subsystems(genome_ids, query_dict, output_file, output_dir, genome_data, s)
    #subsystems_success = run_subsystems_v2(genome_ids, query_dict, output_file, output_dir, genome_data, s)
    proteinfams_success = run_families_v2(genome_ids, query_dict, output_file, output_dir, genome_data, s)

    #generate_report(genome_ids,pathway_success,subsystems_success,proteinfams_success,output_dir)
