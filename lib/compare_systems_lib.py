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

from bvbrc_api import (
    authenticateByEnv,
    getGenomeIdsByGenomeGroup,
    getFeatureDataFrame,
    getSubsystemsDataFrame,
    getPathwayDataFrame,
    getDataForGenomes,
    getQueryData,
    getQueryDataText,
)

import time
import io


def chunker(seq, size):
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def build_key(component_list):
    mod_list = []
    for x in component_list:
        if not isinstance(x, str):
            x = "nan"
        mod_list.append(x)
    # key_parts = [x.lower() for x in mod_list]
    return ":".join(mod_list)


# Returns the entry with the maximum occurences in a pandas dataframe column
# df: pandas dataframe
# col: column name
def get_maximum_value(df, col):
    # value counts returns an ordered list, descending
    col_counts = df[col].value_counts()
    max_label = col_counts[0].index
    return max_label


def return_columns_to_remove(system, columns):
    if system == "subsystems_genes":
        drop_subsystem_columns = [
            "date_inserted",
            "date_modified",
            "genome_name",
            "gene",
            "owner",
            "patric_id",
            "public",
            "product",
            "refseq_locus_tag",
            "taxon_id",
            "_version_",
        ]
        table_columns = list(
            set.intersection(set(drop_subsystem_columns), set(columns))
        )
        return table_columns
    elif system == "subsystems_subsystems":
        drop_subsystem_columns = [
            "feature_id",
            "public",
            "role_id",
            "genome_id",
            "taxon_id",
            "role_name",
            "owner",
            "product",
            "patric_id",
            "genome_name",
            "id",
            "_version_",
            "date_inserted",
            "date_modified",
        ]
        table_columns = list(
            set.intersection(set(drop_subsystem_columns), set(columns))
        )
        return table_columns
    elif system == "proteinfamilies_plfams":
        drop_plfams_columns = [
            "genome_name",
            "accession",
            "patric_id",
            "refseq_locus_tag",
            "alt_locus_tag",
            "feature_id",
            "annotation",
            "feature_type",
            "start",
            "end",
            "strand",
            "figfam_id",
            "pgfam_id",
            "protein_id",
            "aa_length",
            "gene",
            "go",
        ]
        table_columns = list(set.intersection(set(drop_plfams_columns), set(columns)))
        return table_columns
    elif system == "proteinfamilies_pgfams":
        drop_pgfams_columns = [
            "genome_name",
            "accession",
            "patric_id",
            "refseq_locus_tag",
            "alt_locus_tag",
            "feature_id",
            "annotation",
            "feature_type",
            "start",
            "end",
            "strand",
            "figfam_id",
            "plfam_id",
            "protein_id",
            "aa_length",
            "gene",
            "go",
        ]
        table_columns = list(set.intersection(set(drop_pgfams_columns), set(columns)))
        return table_columns
    elif system == "pathways_genes":
        drop_gene_columns = [
            "genome_name",
            "accession",
            "alt_locus_tag",
            "refseq_locus_tag",
            "feature_id",
            "annotation",
            "product",
        ]
        table_columns = list(set.intersection(set(drop_gene_columns), set(columns)))
        return table_columns
    else:  # pathways does not have drop columns
        sys.stderr.write("Error, system is not a valid type\n")
        return []


def run_families(
    genome_ids,
    query_dict,
    output_file,
    output_dir,
    genome_data,
    genome_group_dict,
    session,
):
    print("starting protein families")
    data_dict = {}
    data_dict["plfam"] = {}
    data_dict["pgfam"] = {}
    plfam_genomes = {}
    pgfam_genomes = {}
    present_genome_ids = set()
    genomes_missing_data = {}
    for gids in chunker(genome_ids, 20):
        base = "https://www.bv-brc.org/api/genome_feature/?http_download=true"
        query = f"in(genome_id,({','.join(gids)}))&limit(2500000)&sort(+feature_id)&eq(annotation,PATRIC)"
        headers = {
            "accept": "text/tsv",
            "content-type": "application/rqlquery+x-www-form-urlencoded",
            "Authorization": session.headers["Authorization"],
        }
        result_header = True
        for line in getQueryData(base, query, headers):
            if result_header:
                result_header = False
                print(line)
                continue
            line = line.strip().split("\t")
            # 20 entries in query result with pgfam and plfam data
            if len(line) < 20:
                continue
            try:
                genome_id = line[1].replace('"', "")
                plfam_id = line[14].replace('"', "")
                pgfam_id = line[15].replace('"', "")
                aa_length = line[17].replace('"', "")
                product = line[19].replace('"', "")
            except Exception as e:
                sys.stderr.write(f"Error with the following line:\n{e}\n{line}\n")
                continue
            if aa_length == "":
                continue
            present_genome_ids.add(genome_id)
            ### add to missing genomes data dict
            if genome_id not in genomes_missing_data:
                genomes_missing_data[genome_id] = True
            ### plfam counts
            if plfam_id != "":
                if plfam_id not in data_dict["plfam"]:
                    data_dict["plfam"][plfam_id] = {}
                    data_dict["plfam"][plfam_id]["aa_length_list"] = []
                    data_dict["plfam"][plfam_id]["feature_count"] = 0
                    data_dict["plfam"][plfam_id]["genome_count"] = 0
                    data_dict["plfam"][plfam_id]["product"] = product
                if plfam_id not in plfam_genomes:
                    plfam_genomes[plfam_id] = {}
                if genome_id not in plfam_genomes[plfam_id]:
                    plfam_genomes[plfam_id][genome_id] = 0
                data_dict["plfam"][plfam_id]["aa_length_list"].append(int(aa_length))
                data_dict["plfam"][plfam_id]["feature_count"] += 1
                data_dict["plfam"][plfam_id]["genome_count"] = len(
                    plfam_genomes[plfam_id]
                )
                plfam_genomes[plfam_id][genome_id] += 1
            ### pgfam counts
            if pgfam_id != "":
                if pgfam_id not in data_dict["pgfam"]:
                    data_dict["pgfam"][pgfam_id] = {}
                    data_dict["pgfam"][pgfam_id]["aa_length_list"] = []
                    data_dict["pgfam"][pgfam_id]["feature_count"] = 0
                    data_dict["pgfam"][pgfam_id]["genome_count"] = 0
                    data_dict["pgfam"][pgfam_id]["product"] = product
                if pgfam_id not in pgfam_genomes:
                    pgfam_genomes[pgfam_id] = {}
                if genome_id not in pgfam_genomes[pgfam_id]:
                    pgfam_genomes[pgfam_id][genome_id] = 0
                data_dict["pgfam"][pgfam_id]["aa_length_list"].append(int(aa_length))
                data_dict["pgfam"][pgfam_id]["feature_count"] += 1
                data_dict["pgfam"][pgfam_id]["genome_count"] = len(
                    pgfam_genomes[pgfam_id]
                )
                pgfam_genomes[pgfam_id][genome_id] += 1

    # - get protein family description data
    product_dict = {}
    for plids_list in chunker(list(data_dict["plfam"].keys()), 5000):
        print(f"plids_list has {len(plids_list)} elements")
        base = "https://www.bv-brc.org/api/protein_family_ref/?http_download=true"
        query = (
            f"in(family_id,({','.join(plids_list)}))&limit(2500000)&sort(+family_id)"
        )
        headers = {
            "accept": "application/json",
            "content-type": "application/rqlquery+x-www-form-urlencoded",
            "Authorization": session.headers["Authorization"],
        }
        # headers = {"accept":"text/tsv", "content-type":"application/rqlquery+x-www-form-urlencoded", 'Authorization': session.headers['Authorization']}
        res_data = getQueryDataText(base, query, headers, print_query=False)
        text_data = json.loads(res_data)
        print(f"text_data has {len(text_data)} elements")
        for entry in text_data:
            product_dict[entry["family_id"]] = entry["family_product"]
        if len(text_data) == 0:
            print(query)
    for pgids_list in chunker(list(data_dict["pgfam"].keys()), 5000):
        print(f"pgids_list has {len(pgids_list)} elements")
        base = "https://www.bv-brc.org/api/protein_family_ref/?http_download=true"
        query = (
            f"in(family_id,({','.join(pgids_list)}))&limit(2500000)&sort(+family_id)"
        )
        headers = {
            "accept": "application/json",
            "content-type": "application/rqlquery+x-www-form-urlencoded",
            "Authorization": session.headers["Authorization"],
        }
        res_data = getQueryDataText(base, query, headers, print_query=False)
        text_data = json.loads(res_data)
        print(f"text_data has {len(text_data)} elements")
        for entry in text_data:
            product_dict[entry["family_id"]] = entry["family_product"]
        if len(text_data) == 0:
            print(query)

    # go back and get the mean, max, min, std dev for each family_id
    plfam_line_list = []
    pgfam_line_list = []
    plfam_genome_list = {}
    pgfam_genome_list = {}
    genome_str_dict = {}
    genome_str_dict["plfam"] = {}
    genome_str_dict["pgfam"] = {}
    for plfam_id in data_dict["plfam"]:
        if plfam_id == "":
            continue
        aa_length_list = data_dict["plfam"][plfam_id]["aa_length_list"]
        aa_length_max = max(aa_length_list)
        aa_length_min = min(aa_length_list)
        aa_length_mean = np.mean(aa_length_list)
        aa_length_std = np.std(aa_length_list)
        feature_count = data_dict["plfam"][plfam_id]["feature_count"]
        genome_count = data_dict["plfam"][plfam_id]["genome_count"]
        # genomes = format(feature_count,'#04x').replace('0x','')
        genomes_dir = {}
        plfam_genome_list[plfam_id] = []
        for gid in genome_ids:
            if gid in plfam_genomes[plfam_id]:
                genomes_dir[gid] = format(plfam_genomes[plfam_id][gid], "#04x").replace(
                    "0x", ""
                )
                plfam_genome_list[plfam_id].append(gid)
            else:
                genomes_dir[gid] = "00"
        # genomes = ''.join(genomes_list)
        genome_str_dict["plfam"][plfam_id] = genomes_dir
        # product = data_dict['plfam'][plfam_id]['product']
        if plfam_id in product_dict:
            product = product_dict[plfam_id]
        else:
            product = "NOTHING"
        plfam_str = f"{plfam_id}\t{feature_count}\t{genome_count}\t{product}\t{aa_length_min}\t{aa_length_max}\t{aa_length_mean}\t{aa_length_std}"
        plfam_line_list.append(plfam_str)
    for pgfam_id in data_dict["pgfam"]:
        if pgfam_id == "":
            continue
        aa_length_list = data_dict["pgfam"][pgfam_id]["aa_length_list"]
        aa_length_max = max(aa_length_list)
        aa_length_min = min(aa_length_list)
        aa_length_mean = np.mean(aa_length_list)
        aa_length_std = np.std(aa_length_list)
        feature_count = data_dict["pgfam"][pgfam_id]["feature_count"]
        genome_count = data_dict["pgfam"][pgfam_id]["genome_count"]
        # genomes = format(feature_count,'#04x').replace('0x','')
        genomes_dir = {}
        pgfam_genome_list[pgfam_id] = []
        for gid in genome_ids:
            if gid in pgfam_genomes[pgfam_id]:
                # genomes+=format(pgfam_genomes[pgfam_id][gid],'#04x').replace('0x','')
                genomes_dir[gid] = format(pgfam_genomes[pgfam_id][gid], "#04x").replace(
                    "0x", ""
                )
                pgfam_genome_list[pgfam_id].append(gid)
            else:
                genomes_dir[gid] = "00"
        # genomes = ''.join(genomes_list)
        genome_str_dict["pgfam"][pgfam_id] = genomes_dir
        # product = data_dict['pgfam'][pgfam_id]['product']
        if pgfam_id in product_dict:
            product = product_dict[pgfam_id]
        else:
            product = "NOTHING"
        pgfam_str = f"{pgfam_id}\t{feature_count}\t{genome_count}\t{product}\t{aa_length_min}\t{aa_length_max}\t{aa_length_mean}\t{aa_length_std}"
        pgfam_line_list.append(pgfam_str)

    # output_json['genome_ids'] = genome_ids
    # output_json['genome_ids'] = list(set(genome_ids).intersection(present_genome_ids))

    unsorted_genome_ids = [gid for gid in genome_ids if gid in present_genome_ids]
    # genome_data columns come from the BV-BRC API as lowercase field names
    gid_col = "genome_id" if "genome_id" in genome_data.columns else "Genome ID"
    gname_col = "genome_name" if "genome_name" in genome_data.columns else "Genome Name"
    # Ensure genome ID values are strings — IDs like 292462.10 may be read as
    # float 292462.1 when the dtype hint doesn't match the TSV column name.
    genome_data[gid_col] = genome_data[gid_col].astype(str)
    tmp_data = genome_data.loc[genome_data[gid_col].isin(unsorted_genome_ids)]
    tmp_data = tmp_data.set_index(gid_col)
    tmp_data = tmp_data.loc[unsorted_genome_ids]
    unsorted_genome_names = tmp_data[gname_col].tolist()
    sorted_genome_names, sorted_genome_ids = zip(
        *sorted(zip(unsorted_genome_names, unsorted_genome_ids))
    )

    # add genomes string to each line
    for x in range(0, len(plfam_line_list)):
        line_parts = plfam_line_list[x].split("\t")
        genomes_dir = genome_str_dict["plfam"][line_parts[0]]
        genome_str = ""
        for gid in sorted_genome_ids:
            genome_str += genomes_dir[gid]
        line_parts.append(genome_str)
        plfam_line_list[x] = "\t".join(line_parts)
    for x in range(0, len(pgfam_line_list)):
        line_parts = pgfam_line_list[x].split("\t")
        genomes_dir = genome_str_dict["pgfam"][line_parts[0]]
        genome_str = ""
        for gid in sorted_genome_ids:
            genome_str += genomes_dir[gid]
        line_parts.append(genome_str)
        pgfam_line_list[x] = "\t".join(line_parts)

    header = "family_id\tfeature_count\tgenome_count\tproduct\taa_length_min\taa_length_max\taa_length_mean\taa_length_std\tgenomes"
    plfam_line_list.insert(0, header)
    pgfam_line_list.insert(0, header)

    output_json = {}
    output_json["plfam"] = "\n".join(plfam_line_list)
    output_json["pgfam"] = "\n".join(pgfam_line_list)
    output_json["genome_ids"] = sorted_genome_ids
    output_json["genome_names"] = sorted_genome_names
    output_json["job_name"] = output_file
    output_json["plfam_genomes"] = plfam_genome_list
    output_json["pgfam_genomes"] = pgfam_genome_list

    # add genome groups and other metadata for genome ids to output json
    output_json["genome_data"] = {}
    # Map from output key -> possible column names (API lowercase or display-style)
    extra_fields_map = {
        "isolation_country": ["isolation_country", "Isolation Country"],
        "collection_year": ["collection_year", "Collection Year"],
        "geographic_group": ["geographic_group", "Geographic Group"],
        "host_group": ["host_group", "Host Group"],
        "genome_status": ["genome_status", "Genome Status"],
    }
    for field_key in extra_fields_map:
        output_json["genome_data"][field_key] = []
    gid_col = "genome_id" if "genome_id" in genome_data.columns else "Genome ID"
    tmp_data = genome_data.set_index(gid_col)
    out_genome_groups = []
    for gi in sorted_genome_ids:
        out_genome_groups.append(genome_group_dict[gi])
        for field_key, candidates in extra_fields_map.items():
            # Find the first matching column name
            col = next((c for c in candidates if c in tmp_data.columns), None)
            if col is not None:
                value = tmp_data.loc[gi][col]
            else:
                value = "n/a"
            if isinstance(value, (int, float)) and np.isnan(value):
                value = "n/a"
            elif not isinstance(value, str):
                value = str(value)
            output_json["genome_data"][field_key].append(value)
    output_json["genome_data"]["genome_group"] = out_genome_groups

    output_json_file = os.path.join(
        output_dir, output_file + "_proteinfams_tables.json"
    )
    with open(output_json_file, "w") as o:
        o.write(json.dumps(output_json))

    print("ProteinFamilies Complete")
    return {"success": True, "genomes": present_genome_ids}


def run_subsystems(
    genome_ids, query_dict, output_file, output_dir, genome_data, session
):
    print("starting subsystems")
    subsystems_file = os.path.join(output_dir, output_file + "_subsystems.tsv")
    subsystem_line_list = []
    subsystem_header = (
        "superclass\tclass\tsubclass\tsubsystem_name\tgene_count\trole_count"
    )
    subsystem_line_list.append(subsystem_header)
    subsystem_dict = {}
    overview_counts_dict = {}
    unique_subsystem_features = {}

    subsystem_query_data = []
    required_fields = [
        "superclass",
        "class",
        "subclass",
        "subsystem_name",
        "subsystem_id",
        "feature_id",
        "gene",
        "product",
        "role_id",
        "role_name",
    ]
    subsystem_data_found = False
    subsystem_genomes_found = set()
    subsystem_table_header = None
    print_one = True
    genome_dict = {}
    variant_counts_dict = {}
    genome_data_dict = {}
    for gids in chunker(genome_ids, 20):
        base = "https://www.bv-brc.org/api/subsystem/?http_download=true"
        query = f"in(genome_id,({','.join(gids)}))&limit(2500000)&sort(+id)"
        headers = {
            "accept": "application/json",
            "content-type": "application/rqlquery+x-www-form-urlencoded",
            "Authorization": session.headers["Authorization"],
        }

        # dict_keys(['active', 'class', 'date_inserted', 'date_modified', 'feature_id', 'gene', 'genome_id', 'genome_name', 'id', 'owner', 'patric_id', 'product', 'public', 'refseq_locus_tag', 'role_id', 'role_name', 'subclass', 'subsystem_id', 'subsystem_name', 'superclass', 'taxon_id', '_version_'])

        print("Query = {0}\nHeaders = {1}".format(base + "&" + query, headers))
        result_header = True
        current_header = None
        all_data = json.loads(getQueryDataText(base, query, headers))
        for line in all_data:
            subsystem_data_found = True
            if result_header:
                result_header = False
                print(line)
                current_header = line.keys()
                if subsystem_table_header is None or len(current_header) > len(
                    subsystem_table_header
                ):
                    # subsystem_table_header = current_header
                    subsystem_table_header = [x for x in current_header]
                continue
            if print_one:
                print_one = False
                print(line)
            subsystem_fields = line
            for field in required_fields:
                if field not in subsystem_fields:
                    subsystem_fields[field] = ""
            subsystem_query_data.append(subsystem_fields)
            try:
                active = subsystem_fields["active"]
                clss = subsystem_fields["class"]
                feature_id = subsystem_fields["feature_id"]
                gene = subsystem_fields["gene"]
                genome_id = subsystem_fields["genome_id"]
                genome_name = subsystem_fields["genome_name"]
                product = subsystem_fields["product"]
                role_id = subsystem_fields["role_id"]
                role_name = subsystem_fields["role_name"]
                subclass = subsystem_fields["subclass"]
                subsystem_id = subsystem_fields["subsystem_id"]
                subsystem_name = subsystem_fields["subsystem_name"]
                superclass = subsystem_fields["superclass"]
                # TODO: underlying issue of metadata, capitalizition of superclasses is not consistent
                superclass = superclass.upper()
            except Exception as e:
                sys.stderr.write(f"Error with the following line:\n{e}\n{line}\n")
                continue
            subsystem_genomes_found.add(genome_id)
            if genome_id not in genome_dict:
                genome_dict[genome_name] = genome_id
                genome_data_dict[genome_id] = {}
                genome_data_dict[genome_id]["genes"] = []
            if superclass not in subsystem_dict:
                subsystem_dict[superclass] = {}
                overview_counts_dict[superclass] = {}
            if clss not in subsystem_dict[superclass]:
                subsystem_dict[superclass][clss] = {}
                overview_counts_dict[superclass][clss] = {}
            if subclass not in subsystem_dict[superclass][clss]:
                subsystem_dict[superclass][clss][subclass] = {}
                overview_counts_dict[superclass][clss][subclass] = {}
                overview_counts_dict[superclass][clss][subclass]["subsystem_names"] = (
                    set()
                )
                overview_counts_dict[superclass][clss][subclass]["gene_set"] = set()
            if subsystem_name not in subsystem_dict[superclass][clss][subclass]:
                subsystem_dict[superclass][clss][subclass][subsystem_name] = {}
                subsystem_dict[superclass][clss][subclass][subsystem_name][
                    "gene_set"
                ] = set()
                subsystem_dict[superclass][clss][subclass][subsystem_name][
                    "role_set"
                ] = set()
                subsystem_dict[superclass][clss][subclass][subsystem_name][
                    "active_genome_dict"
                ] = {}
                subsystem_dict[superclass][clss][subclass][subsystem_name][
                    "subsystem_id"
                ] = subsystem_id
                subsystem_dict[superclass][clss][subclass][subsystem_name][
                    "subsystem_counts"
                ] = 0
            overview_counts_dict[superclass][clss][subclass]["subsystem_names"].add(
                subsystem_name
            )
            subsystem_dict[superclass][clss][subclass][subsystem_name][
                "active_genome_dict"
            ][genome_id] = active
            # sub_key = superclass + clss + subclass + subsystem_name
            if subsystem_id not in variant_counts_dict:
                variant_counts_dict[subsystem_id] = {}
                variant_counts_dict[subsystem_id]["active"] = 0
                variant_counts_dict[subsystem_id]["likely"] = 0
                variant_counts_dict[subsystem_id]["inactive"] = 0
            if active == "active" or active == "likely":
                variant_counts_dict[subsystem_id][active] += 1
            else:  # never reached, the genome just doesn't have an entry
                variant_counts_dict[subsystem_id]["inactive"] += 1
            # TODO: repeated features; do I count these towards the gene counts?
            # if feature_id in subsystem_dict[superclass][clss][subclass][subsystem_name]['gene_set']:
            #    with open('repeated_feature_ids.txt','a') as o:
            #       o.write(f'{feature_id}\n')
            # if feature_id not in genome_data_dict[genome_id]["genes"]:
            # subsystem_dict[superclass][clss][subclass][subsystem_name]['gene_set'].add(feature_id)
            subsystem_dict[superclass][clss][subclass][subsystem_name]["gene_set"].add(
                feature_id
            )
            # genome_data_dict[genome_id]["genes"].append(feature_id)
            overview_counts_dict[superclass][clss][subclass]["gene_set"].add(feature_id)
            # if role_id is not None and role_id != '':
            subsystem_dict[superclass][clss][subclass][subsystem_name]["role_set"].add(
                role_id
            )
            subsystem_dict[superclass][clss][subclass][subsystem_name][
                "subsystem_counts"
            ] += 1

    if not subsystem_data_found:
        return {"success": False}

    parsed_query_data = []
    for line in subsystem_query_data:
        new_line = ""
        for field in subsystem_table_header:
            if new_line != "":
                new_line += "\t"
            if field not in line:
                new_line += " "
            else:
                value = line[field]
                if not isinstance(value, str):
                    value = str(value)
                new_line += value
        parsed_query_data.append(new_line.split("\t"))

    subsystem_df = pd.DataFrame(parsed_query_data, columns=subsystem_table_header)
    subsystem_df.to_csv(subsystems_file, index=False, sep="\t")

    gene_df = query_dict["feature"]

    # Ensure merge keys are string type — genome_id can be inferred as
    # float64 by pd.read_csv when all values in a batch look numeric.
    for merge_col in ["genome_id", "feature_id"]:
        if merge_col in gene_df.columns:
            gene_df[merge_col] = gene_df[merge_col].astype(str)
        if merge_col in subsystem_df.columns:
            subsystem_df[merge_col] = subsystem_df[merge_col].astype(str)

    gene_df = pd.merge(
        gene_df,
        subsystem_df.drop(
            return_columns_to_remove("subsystems_genes", subsystem_df.columns.tolist()),
            axis=1,
        ),
        on=["genome_id", "feature_id"],
        how="inner",
    )

    # get data for conservation scores
    unique_subsystem_features = {}
    unique_subsystem_roles = {}
    for idx in range(0, gene_df.shape[0]):
        curr_sub_id = gene_df.iloc[idx].subsystem_id
        if curr_sub_id not in unique_subsystem_features:
            unique_subsystem_features[curr_sub_id] = {}
        if curr_sub_id not in unique_subsystem_roles:
            unique_subsystem_roles[curr_sub_id] = {}
        gene = gene_df.iloc[idx]["gene"]
        role = gene_df.iloc[idx]["role_id"]
        genome_id = gene_df.iloc[idx].genome_id
        if not gene is None and not gene is np.nan:
            if gene not in unique_subsystem_features[curr_sub_id]:
                unique_subsystem_features[curr_sub_id][gene] = set()
            unique_subsystem_features[curr_sub_id][gene].add(genome_id)
        if not role is None and not role is np.nan:
            if role not in unique_subsystem_roles[curr_sub_id]:
                unique_subsystem_roles[curr_sub_id][role] = set()
            unique_subsystem_roles[curr_sub_id][role].add(genome_id)

    # conservation scores
    # gets counts for overview dict, any other adjustments
    # create subsystems table
    subsystems_table_list = []
    overview_dict = {}
    for superclass in overview_counts_dict:
        overview_dict[superclass] = {}
        overview_dict[superclass]["subsystem_name_counts"] = 0
        overview_dict[superclass]["gene_counts"] = 0
        for clss in overview_counts_dict[superclass]:
            overview_dict[superclass][clss] = {}
            overview_dict[superclass][clss]["subsystem_name_counts"] = 0
            overview_dict[superclass][clss]["gene_counts"] = 0
            for subclass in overview_counts_dict[superclass][clss]:
                overview_dict[superclass][clss][subclass] = {}
                overview_dict[superclass][clss][subclass]["subsystem_name_counts"] = (
                    len(
                        overview_counts_dict[superclass][clss][subclass][
                            "subsystem_names"
                        ]
                    )
                )
                overview_dict[superclass][clss][subclass]["gene_counts"] = len(
                    overview_counts_dict[superclass][clss][subclass]["gene_set"]
                )
                overview_dict[superclass][clss]["subsystem_name_counts"] += len(
                    overview_counts_dict[superclass][clss][subclass]["subsystem_names"]
                )
                overview_dict[superclass][clss]["gene_counts"] += len(
                    overview_counts_dict[superclass][clss][subclass]["gene_set"]
                )
                overview_dict[superclass]["gene_counts"] += len(
                    overview_counts_dict[superclass][clss][subclass]["gene_set"]
                )
                overview_dict[superclass]["subsystem_name_counts"] += len(
                    overview_counts_dict[superclass][clss][subclass]["subsystem_names"]
                )
                for subsystem_name in subsystem_dict[superclass][clss][subclass]:
                    # sub_key = superclass + clss + subclass + subsystem_name
                    subsystem_id = subsystem_dict[superclass][clss][subclass][
                        subsystem_name
                    ]["subsystem_id"]
                    # gene conservation
                    # if subsystem_id in unique_subsystem_features:
                    # for gene in unique_subsystem_features[subsystem_id]:
                    #    gene_numerator += len(unique_subsystem_features[subsystem_id][gene])
                    gene_numerator = len(
                        subsystem_dict[superclass][clss][subclass][subsystem_name][
                            "gene_set"
                        ]
                    )
                    gene_denominator = len(
                        subsystem_dict[superclass][clss][subclass][subsystem_name][
                            "role_set"
                        ]
                    ) * len(subsystem_genomes_found)
                    gene_conservation = 0
                    if gene_denominator > 0:
                        gene_conservation = float(gene_numerator) / float(
                            gene_denominator
                        )
                    # role conservation
                    role_numerator = 0
                    role_denominator = 0
                    if subsystem_id in unique_subsystem_roles:
                        for role in unique_subsystem_roles[subsystem_id]:
                            role_numerator += len(
                                unique_subsystem_roles[subsystem_id][role]
                            )
                        role_denominator = len(
                            subsystem_dict[superclass][clss][subclass][subsystem_name][
                                "role_set"
                            ]
                        ) * len(subsystem_genomes_found)
                    role_conservation = 0
                    if role_denominator > 0:
                        role_conservation = (
                            float(role_numerator) / float(role_denominator) * 100
                        )
                    active_num = len(
                        subsystem_df[
                            subsystem_df.subsystem_id == subsystem_id
                        ].genome_id.unique()
                    )
                    new_entry = {
                        "superclass": superclass,
                        "class": clss,
                        "subclass": subclass,
                        "subsystem_name": subsystem_name,
                        "subsystem_id": subsystem_id,
                        "role_counts": len(
                            subsystem_dict[superclass][clss][subclass][subsystem_name][
                                "role_set"
                            ]
                        ),
                        "gene_counts": len(
                            subsystem_dict[superclass][clss][subclass][subsystem_name][
                                "gene_set"
                            ]
                        ),
                        "genome_count": len(
                            subsystem_dict[superclass][clss][subclass][subsystem_name][
                                "active_genome_dict"
                            ]
                        ),
                        "gene_conservation": gene_conservation,
                        "role_conservation": role_conservation,
                        "prop_active": float(active_num)
                        / float(len(subsystem_genomes_found)),
                    }
                    subsystems_table_list.append(new_entry)
    subsystems_table = pd.DataFrame(subsystems_table_list)

    # Variant matrix
    variant_mtx_header = "\t\t\t\t\t\t"
    gid_str = ""
    genome_name_list = list(genome_dict.keys())
    genome_name_list.sort()
    for genome_name in genome_name_list:
        variant_mtx_header += f"\t{genome_name}"
        gid_str += f"\t{genome_dict[genome_name]}"
    variant_mtx_header += "\nSuperclass\tClass\tSubclass\tSS\tactive\tlikely\tinactive"
    variant_mtx_header += gid_str
    variant_mtx_lines = []
    variant_mtx_lines.append(variant_mtx_header)
    for superclass in subsystem_dict:
        for clss in subsystem_dict[superclass]:
            for subclass in subsystem_dict[superclass][clss]:
                for subsystem_name in subsystem_dict[superclass][clss][subclass]:
                    # sub_key = superclass + clss + subclass + subsystem_name
                    subsystem_id = subsystem_dict[superclass][clss][subclass][
                        subsystem_name
                    ]["subsystem_id"]
                    inactive_value = 0
                    new_var_line_p1 = (
                        f"{superclass}\t{clss}\t{subclass}\t{subsystem_name}"
                    )
                    new_var_line_p1 += f"\t{variant_counts_dict[subsystem_id]['active']}\t{variant_counts_dict[subsystem_id]['likely']}\t{inactive_value}"
                    new_var_line_p2 = ""
                    for genome_name in genome_name_list:
                        if (
                            genome_dict[genome_name]
                            in subsystem_dict[superclass][clss][subclass][
                                subsystem_name
                            ]["active_genome_dict"]
                        ):
                            new_var_line_p2 += f"\t{subsystem_dict[superclass][clss][subclass][subsystem_name]['active_genome_dict'][genome_dict[genome_name]]}"
                        else:
                            new_var_line_p2 += f"\tinactive"
                            inactive_value += 1
                    new_var_line = (
                        new_var_line_p1 + f"\t{inactive_value}" + new_var_line_p2
                    )
                    variant_mtx_lines.append(new_var_line)
    variant_mtx_text = "\n".join(variant_mtx_lines)
    variant_mtx_file = subsystems_file.replace(".tsv", "_variant_mtx.tsv")
    with open(variant_mtx_file, "w") as o:
        o.write(variant_mtx_text)

    output_json_file = subsystems_file.replace(".tsv", "_tables.json")

    output_json = {}
    output_json["genome_ids"] = list(subsystem_genomes_found)
    # genome_data columns come from the BV-BRC API as lowercase field names
    gid_col = "genome_id" if "genome_id" in genome_data.columns else "Genome ID"
    gname_col = "genome_name" if "genome_name" in genome_data.columns else "Genome Name"
    # Ensure genome ID values are strings — the API returns TSV with display
    # headers ("Genome ID" not "genome_id") so the dtype hint in read_csv has
    # no effect and IDs like 292462.10 get read as float 292462.1.
    genome_data[gid_col] = genome_data[gid_col].astype(str)
    output_json["genome_names"] = (
        genome_data.set_index(gid_col)
        .loc[list(subsystem_genomes_found)][gname_col]
        .tolist()
    )  # returns a list of genome names in the same order as the genome ids
    output_json["overview"] = overview_dict
    output_json["job_name"] = output_file
    output_json["subsystems"] = subsystems_table.to_csv(index=False, sep="\t")
    output_json["genes"] = gene_df.to_csv(index=False, sep="\t")
    with open(output_json_file, "w") as o:
        o.write(json.dumps(output_json))

    print("Subsystems complete")
    return {"success": True, "genomes": list(subsystem_genomes_found)}


def run_pathways(genome_ids, query_dict, output_file, output_dir, genome_data, session):
    print("starting pathways")
    pathways_file = os.path.join(output_dir, output_file + "_pathways.tsv")
    # pathway_df = query_dict['pathway']
    # pathway_df.to_csv(pathways_file,sep='\t',index=False)
    pathway_line_list = []
    ec_line_list = []
    pathway_header = "annotation\tpathway_id\tpathway_name\tpathway_class\tgenome_count\tec_count\tgene_count\tgenome_ec\tec_conservation\tgene_conservation"
    ec_header = "annotation\tpathway_id\tpathway_name\tpathway_class\tec_description\tec_number\tgenome_count\tec_count\tgene_count\tgenome_ec"
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
    required_fields = [
        "annotation",
        "ec_description",
        "ec_number",
        "feature_id",
        "genome_id",
        "pathway_class",
        "pathway_id",
        "pathway_name",
        "patric_id",
        "product",
    ]
    pathway_data_found = False
    pathway_genomes_found = set()
    pathway_table_header = None
    for gids in chunker(genome_ids, 20):
        base = "https://www.bv-brc.org/api/pathway/?http_download=true"
        query = f"in(genome_id,({','.join(gids)}))&limit(2500000)&sort(+id)&eq(annotation,PATRIC)"
        headers = {
            "accept": "application/json",
            "content-type": "application/rqlquery+x-www-form-urlencoded",
            "Authorization": session.headers["Authorization"],
        }

        print("Query = {0}\nHeaders = {1}".format(base + "&" + query, headers))
        # accession       alt_locus_tag   annotation      date_inserted   date_modified   ec_description  ec_number       feature_id      genome_ec       genome_id       genome_name     id      owner   pathway_class   pathway_ec      pathway_id   pathway_name     patric_id       product public  refseq_locus_tag        sequence_id     taxon_id        _version_

        result_header = True
        current_header = None
        all_data = json.loads(getQueryDataText(base, query, headers))
        for line in all_data:
            pathway_data_found = True
            if result_header:
                result_header = False
                print(line)
                current_header = line.keys()
                if pathway_table_header is None or len(current_header) > len(
                    pathway_table_header
                ):
                    pathway_table_header = current_header
                # pathway_query_data.append(line)
                continue
            # line = line.strip().replace('\"','').split('\t')
            pathway_fields = line
            # for idx,f in enumerate(line):
            #    pathway_fields[pathway_table_header[idx]] = f
            for field in required_fields:
                if field not in pathway_fields:
                    pathway_fields[field] = ""
            pathway_query_data.append(pathway_fields)
            try:
                annotation = pathway_fields["annotation"]
                ec_description = pathway_fields["ec_description"]
                ec_number = pathway_fields["ec_number"]
                feature_id = pathway_fields["feature_id"]
                genome_id = pathway_fields["genome_id"]
                # genome_name = line[10].replace('\"','')
                pathway_class = pathway_fields["pathway_class"]
                pathway_id = pathway_fields["pathway_id"]
                pathway_name = pathway_fields["pathway_name"]
                patric_id = pathway_fields["patric_id"]
                product = pathway_fields["product"]
            except Exception as e:
                sys.stderr.write(f"Error with the following line:\n{e}\n{line}\n")
                continue
            pathway_genomes_found.add(genome_id)
            unique_pathways.add(pathway_id)
            unique_ecs.add(ec_number)
            # unique_features.add(patric_id)
            """
            if pathway_id not in unique_pathway_features:
                unique_pathway_features[pathway_id] = {} 
            if 'gene' in pathway_fields:
                pathway_gene = pathway_fields['gene']
                if pathway_gene not in unique_pathway_features[pathway_id]: 
                    unique_pathway_features[pathway_id][pathway_gene] = set()
                unique_pathway_features[pathway_id][pathway_gene].add(genome_id)
            """
            if pathway_id not in unique_pathway_ecs:
                unique_pathway_ecs[pathway_id] = {}
            if ec_number not in unique_pathway_ecs[pathway_id]:
                unique_pathway_ecs[pathway_id][ec_number] = set()
            unique_pathway_ecs[pathway_id][ec_number].add(genome_id)

            # pathway data
            if pathway_id not in pathway_dict:
                pathway_dict[pathway_id] = {}
                pathway_dict[pathway_id]["annotation"] = annotation
                pathway_dict[pathway_id]["pathway_id"] = pathway_id
                pathway_dict[pathway_id]["pathway_name"] = pathway_name
                pathway_dict[pathway_id]["pathway_class"] = pathway_class
                pathway_dict[pathway_id]["genome_count"] = set()
                pathway_dict[pathway_id]["ec_count"] = set()
                pathway_dict[pathway_id]["gene_count"] = set()
                pathway_dict[pathway_id]["genome_ec"] = set()
            pathway_dict[pathway_id]["genome_count"].add(genome_id)
            pathway_dict[pathway_id]["ec_count"].add(ec_number)
            pathway_dict[pathway_id]["gene_count"].add(feature_id)
            pathway_dict[pathway_id]["genome_ec"].add(genome_id + "_" + ec_number)
            # ec data
            # ec_header = 'annotation\tpathway_id\tpathway_name\tpathway_class\tproduct\tec_number\tgenome_count\tec_count\tgene_count\tgenome_ec'
            if pathway_id not in ec_dict:
                ec_dict[pathway_id] = {}
            if ec_number not in ec_dict[pathway_id]:
                ec_dict[pathway_id][ec_number] = {}
                ec_dict[pathway_id][ec_number]["annotation"] = annotation
                ec_dict[pathway_id][ec_number]["pathway_id"] = pathway_id
                ec_dict[pathway_id][ec_number]["pathway_name"] = pathway_name
                ec_dict[pathway_id][ec_number]["pathway_class"] = pathway_class
                ec_dict[pathway_id][ec_number]["ec_description"] = ec_description
                ec_dict[pathway_id][ec_number]["ec_number"] = ec_number
                ec_dict[pathway_id][ec_number]["genome_count"] = set()
                ec_dict[pathway_id][ec_number]["ec_count"] = set()
                ec_dict[pathway_id][ec_number]["gene_count"] = set()
                ec_dict[pathway_id][ec_number]["genome_ec"] = set()
            ec_dict[pathway_id][ec_number]["genome_count"].add(genome_id)
            ec_dict[pathway_id][ec_number]["ec_count"].add(ec_number)
            ec_dict[pathway_id][ec_number]["gene_count"].add(feature_id)
            ec_dict[pathway_id][ec_number]["genome_ec"].add(genome_id + "_" + ec_number)

    if not pathway_data_found:
        return {"success": False}

    parsed_query_data = []
    for line in pathway_query_data:
        new_line = ""
        for field in pathway_table_header:
            if new_line != "":
                new_line += "\t"
            if field not in line:
                new_line += " "
            else:
                value = line[field]
                if not isinstance(value, str):
                    value = str(value)
                new_line += value
        parsed_query_data.append(new_line.split("\t"))

    pathway_df = pd.DataFrame(parsed_query_data, columns=pathway_table_header)
    gene_df = query_dict["feature"]

    # Ensure merge keys are string type — genome_id can be inferred as
    # float64 by pd.read_csv when all values in a batch look numeric.
    for merge_col in ["genome_id", "patric_id"]:
        if merge_col in gene_df.columns:
            gene_df[merge_col] = gene_df[merge_col].astype(str)
        if merge_col in pathway_df.columns:
            pathway_df[merge_col] = pathway_df[merge_col].astype(str)

    genes_output = pd.merge(
        gene_df.drop(
            return_columns_to_remove("pathways_genes", gene_df.columns.tolist()), axis=1
        ),
        pathway_df,
        on=["genome_id", "patric_id"],
        how="inner",
    )

    if "gene_x" in genes_output.columns:
        genes_output["gene"] = genes_output["gene_x"]
        genes_output.drop(["gene_x", "gene_y"], inplace=True, axis=1)

    for idx in range(0, genes_output.shape[0]):
        pathway_id = genes_output.iloc[idx].pathway_id
        if pathway_id not in unique_pathway_features:
            unique_pathway_features[pathway_id] = {}
        gene = genes_output.iloc[idx]["gene"]
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
        pathway_dict[pathway_id]["genome_count"] = len(
            pathway_dict[pathway_id]["genome_count"]
        )
        pathway_dict[pathway_id]["ec_count"] = len(pathway_dict[pathway_id]["ec_count"])
        pathway_dict[pathway_id]["gene_count"] = len(
            pathway_dict[pathway_id]["gene_count"]
        )
        pathway_dict[pathway_id]["genome_ec"] = len(
            pathway_dict[pathway_id]["genome_ec"]
        )
        # pathway_dict[pathway_id]['ec_conservation'] = float(len(unique_pathway_ecs[pathway_id]))/float(len(unique_ecs))*100.0
        # pathway_dict[pathway_id]['gene_conservation'] = float(len(unique_pathway_features[pathway_id]))/float(len(unique_features))*100.0
        annotation = pathway_dict[pathway_id]["annotation"]
        pathway_id = pathway_dict[pathway_id]["pathway_id"]
        pathway_name = pathway_dict[pathway_id]["pathway_name"]
        pathway_class = pathway_dict[pathway_id]["pathway_class"]
        genome_count = pathway_dict[pathway_id]["genome_count"]
        ec_count = pathway_dict[pathway_id]["ec_count"]
        gene_count = pathway_dict[pathway_id]["gene_count"]
        genome_ec = pathway_dict[pathway_id]["genome_ec"]
        # calculate ec_conservation score
        ec_numerator = 0
        ec_denominator = 0
        for ec_number in unique_pathway_ecs[pathway_id]:
            ec_numerator += len(unique_pathway_ecs[pathway_id][ec_number])
            ec_denominator += len(pathway_genomes_found)
        # ec_numerator = float(ec_numerator) * float(len(pathway_genomes_found))
        # ec_denominator = float(ec_denominator) * float(len(pathway_genomes_found))
        if ec_denominator == 0:
            ec_conservation = 0
        else:
            ec_conservation = float(ec_numerator) / float(ec_denominator) * 100.0
        # calculate gene_conservation
        gene_numerator = gene_count
        gene_denominator = ec_count * len(pathway_genomes_found)
        if gene_denominator == 0:
            gene_conservation = 0
        else:
            gene_conservation = float(gene_numerator) / float(gene_denominator)
        pathway_line = f"{annotation}\t{pathway_id}\t{pathway_name}\t{pathway_class}\t{genome_count}\t{ec_count}\t{gene_count}\t{genome_ec}\t{ec_conservation}\t{gene_conservation}"
        pathway_line_list.append(pathway_line)
        # now EC data
        for ec_number in ec_dict[pathway_id]:
            ec_dict[pathway_id][ec_number]["genome_count"] = len(
                ec_dict[pathway_id][ec_number]["genome_count"]
            )
            ec_dict[pathway_id][ec_number]["ec_count"] = len(
                ec_dict[pathway_id][ec_number]["ec_count"]
            )
            ec_dict[pathway_id][ec_number]["gene_count"] = len(
                ec_dict[pathway_id][ec_number]["gene_count"]
            )
            ec_dict[pathway_id][ec_number]["genome_ec"] = len(
                ec_dict[pathway_id][ec_number]["genome_ec"]
            )
            annotation = ec_dict[pathway_id][ec_number]["annotation"]
            pathway_id = ec_dict[pathway_id][ec_number]["pathway_id"]
            pathway_name = ec_dict[pathway_id][ec_number]["pathway_name"]
            pathway_class = ec_dict[pathway_id][ec_number]["pathway_class"]
            ec_description = ec_dict[pathway_id][ec_number]["ec_description"]
            ec_number = ec_dict[pathway_id][ec_number]["ec_number"]
            genome_count = ec_dict[pathway_id][ec_number]["genome_count"]
            ec_count = ec_dict[pathway_id][ec_number]["ec_count"]
            gene_count = ec_dict[pathway_id][ec_number]["gene_count"]
            genome_ec = ec_dict[pathway_id][ec_number]["genome_ec"]
            ec_line = f"{annotation}\t{pathway_id}\t{pathway_name}\t{pathway_class}\t{ec_description}\t{ec_number}\t{genome_count}\t{ec_count}\t{gene_count}\t{genome_ec}"
            ec_line_list.append(ec_line)

    pathway_output = "\n".join(pathway_line_list)
    ec_output = "\n".join(ec_line_list)

    output_json = {}
    output_json["pathway"] = pathway_output
    output_json["ecnumber"] = ec_output
    output_json["genes"] = genes_output.to_csv(index=False, sep="\t")
    output_json["genome_ids"] = list(pathway_genomes_found)
    output_json["job_name"] = output_file

    pathway_df.to_csv(pathways_file, sep="\t", index=False)

    output_json_file = pathways_file.replace(".tsv", "_tables.json")
    with open(output_json_file, "w") as o:
        o.write(json.dumps(output_json))

    print("Pathways Complete")
    pathway_success_json = {"genomes": list(pathway_genomes_found), "success": True}
    return pathway_success_json


def generate_report(
    genome_ids, pathway_obj, subsystems_obj, proteinfams_obj, output_dir
):
    report_text_list = []
    if pathway_obj["success"]:
        report_text_list.append(
            f"Pathways succeded: {len(pathway_obj['genomes'])} out of {len(genome_ids)} genomes had pathway data"
        )
        if len(pathway_obj["genomes"]) != len(genome_ids):
            missing_pathway_genomes = list(
                set(genome_ids).difference(set(pathway_obj["genomes"]))
            )
            report_text_list.append(
                f"Genomes Missing from Pathways: {','.join(missing_pathway_genomes)}"
            )
    else:
        report_text_list.append("Pathways Failed: see stdout and stderr")
    if subsystems_obj["success"]:
        report_text_list.append(
            f"Subsystems succeeded: {len(subsystems_obj['genomes'])} out of {len(genome_ids)} genomes had subsystems data"
        )
        if len(subsystems_obj["genomes"]) != len(genome_ids):
            missing_subsystems_genomes = list(
                set(genome_ids).difference(set(subsystems_obj["genomes"]))
            )
            report_text_list.append(
                f"Genomes Missing from Subsystems: {','.join(missing_subsystems_genomes)}"
            )
    else:
        report_text_list.append("Subsystems Failed: see stdout and stderr")
    if proteinfams_obj["success"]:
        report_text_list.append(
            f"ProteinFamilies succeeded: {len(proteinfams_obj['genomes'])} out of {len(genome_ids)} genomes had proteinfamilies data"
        )
        if len(proteinfams_obj["genomes"]) != len(genome_ids):
            missing_proteinfams_genomes = list(
                set(genome_ids).difference(set(proteinfams_obj["genomes"]))
            )
            report_text_list.append(
                f"Genomes Missing from ProteinFamilies: {','.join(missing_proteinfams_genomes)}"
            )
    else:
        report_text_list.append("ProteinFamilies Failed: see stdout and sterr")
    report_text = "\n".join(report_text_list)
    report_file = os.path.join(output_dir, "report.txt")
    with open(report_file, "w") as o:
        o.write(report_text)


# Store pathways, subsystems, and features queries in a dictionary
def run_feature_queries(genome_ids, session):
    query_dict = {}
    ### Run features query
    if True:
        print("features query")
        try:
            feature_df = getFeatureDataFrame(genome_ids, session, limit=2500000)
        except Exception as e:
            print(f"Error running features query:\n{e}\n")
            return None
        if not feature_df is None:
            column_map = {
                "Genome": "genome_name",
                "Genome ID": "genome_id",
                "Accession": "accession",
                "BRC ID": "patric_id",
                "RefSeq Locus Tag": "refseq_locus_tag",
                "Alt Locus Tag": "alt_locus_tag",
                "Feature ID": "feature_id",
                "Annotation": "annotation",
                "Feature Type": "feature_type",
                "Start": "start",
                "End": "end",
                "Length": "length",
                "Strand": "strand",
                "FIGfam ID": "figfam_id",
                "PATRIC genus-specific families (PLfams)": "plfam_id",
                "PATRIC cross-genus families (PGfams)": "pgfam_id",
                "Protein ID": "protein_id",
                "AA Length": "aa_length",
                "Gene Symbol": "gene",
                "Product": "product",
                "GO": "go",
            }
            if "Genome ID" in feature_df.columns:
                feature_df.rename(columns=column_map, inplace=True)
            feature_df["plfam_index"] = feature_df["plfam_id"]
            feature_df["pgfam_index"] = feature_df["pgfam_id"]
            # feature_df.set_index('plfam_index', inplace=True)
            query_dict["feature"] = feature_df
        else:
            sys.stderr.write("Features dataframe is None\n")
    return query_dict


def get_genome_group_ids(group_list, session):
    genome_group_ids = []
    genome_group_list = []
    for genome_group in group_list:
        genome_id_list = getGenomeIdsByGenomeGroup(
            genome_group, session, genomeGroupPath=True
        )
        genome_group_ids = genome_group_ids + genome_id_list
        genome_group_list += [genome_group] * len(genome_id_list)
    return (genome_group_ids, genome_group_list)


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

    # TODO: Testing adding genome groups to genomeData
    genome_ids = job_data["genome_ids"]
    genome_group_list = ["None"] * len(genome_ids)
    if len(job_data["genome_groups"]) > 0:
        genome_group_ids, curr_genome_group_list = get_genome_group_ids(
            job_data["genome_groups"], s
        )
        if len(genome_group_ids) == 0:
            sys.stderr.write("FAILED to get genome ids for genome groups: exiting")
            sys.exit(-1)
        # make ids unique
        genome_ids = genome_ids + genome_group_ids
        genome_group_list += curr_genome_group_list

    # create genome group dictionary
    genome_group_dict = {}
    for idx, gi in enumerate(genome_ids):
        gg = os.path.basename(genome_group_list[idx])
        if gi in genome_group_dict:
            genome_group_dict[gi] = genome_group_dict[gi] + "," + gg
        else:
            genome_group_dict[gi] = gg

    genome_ids = list(set(genome_ids))
    # optionally add more genome info to output
    genome_data = getDataForGenomes(genome_ids, s)

    # Ensure the genome ID column is string type.  The BV-BRC API returns TSV
    # with display headers ("Genome ID") but getDataForGenomes uses
    # dtype={'genome_id':str} which only matches API-style headers.  When the
    # header is "Genome ID", IDs like 292462.10 are read as float 292462.1,
    # breaking all downstream lookups that compare against string genome IDs.
    if genome_data is not None:
        gid_col = "genome_id" if "genome_id" in genome_data.columns else "Genome ID"
        genome_data[gid_col] = genome_data[gid_col].astype(str)

    query_dict = run_feature_queries(genome_ids, s)
    if not query_dict:
        sys.stderr.write("Error running features queries: terminating\n")
        report_text = "Error running features queries: see stdout and stderr"
        report_file = os.path.join(output_dir, "report.txt")
        with open(report_file, "w") as o:
            o.write(report_text)
        sys.exit(0)

    # TODO: add multithreading
    pool = multiprocessing.Pool(processes=3)
    pathway_result = pool.apply_async(
        run_pathways,
        args=(genome_ids, query_dict, output_file, output_dir, genome_data, s),
    )
    subsystems_result = pool.apply_async(
        run_subsystems,
        args=(genome_ids, query_dict, output_file, output_dir, genome_data, s),
    )
    proteinfams_result = pool.apply_async(
        run_families,
        args=(
            genome_ids,
            query_dict,
            output_file,
            output_dir,
            genome_data,
            genome_group_dict,
            s,
        ),
    )

    pathway_success = pathway_result.get()
    subsystems_success = subsystems_result.get()
    proteinfams_success = proteinfams_result.get()
    pool.close()
    pool.join()

    generate_report(
        genome_ids, pathway_success, subsystems_success, proteinfams_success, output_dir
    )
