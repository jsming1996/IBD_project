#!/usr/bin/env python3
import sys
import os
import pandas as pd
from collections import defaultdict

FP_feature_dir = sys.argv[1]
SNP_csv = sys.argv[2]

# identify gene id
def add_gene_id(row, FP_feature_list):
	if row["species"] in FP_feature_list.keys():
		species_df = FP_feature_list[row["species"]]
		if len(species_df) == 0:
			return ""
		temp_species_df = species_df.loc[species_df["start"] <= row["snp"]]
		if len(temp_species_df) == 0:
			return ""
		temp_species_df = temp_species_df.loc[species_df["end"] >= row["snp"]]
		if len(temp_species_df) == 0:
			return ""
		temp_species_df1 = temp_species_df.loc[temp_species_df["type"].str.contains("gene")]
		temp_species_df2 = temp_species_df.loc[temp_species_df["type"].str.contains("pseudogenic_tRNA")]
		if len(temp_species_df1) == 0 and len(temp_species_df2) == 0:
			return ""
		elif len(temp_species_df2) == 0:
			return temp_species_df1.iloc[0]["id"]
		else:
			return temp_species_df2.iloc[0]["id"]
	else:
		return ""

# identify gene start position
def add_start(row, FP_feature_list):
	if row["gene_id"] == "":
		return ""
	else:
		species_df = FP_feature_list[row["species"]]
		temp_species_df = species_df.loc[species_df["id"] == row["gene_id"]]
		return temp_species_df.iloc[0]["start"]

# identify gene end position
def add_end(row, FP_feature_list):
	if row["gene_id"] == "":
		return ""
	else:
		species_df = FP_feature_list[row["species"]]
		temp_species_df = species_df.loc[species_df["id"] == row["gene_id"]]
		return temp_species_df.iloc[0]["end"]

# identify protein id
def add_protein_id(row, FP_feature_list):
	if row["species"] in FP_feature_list.keys():
		species_df = FP_feature_list[row["species"]]
		if len(species_df) == 0:
			return ""
		temp_species_df = species_df.loc[species_df["start"] <= row["snp"]]
		if len(temp_species_df) == 0:
			return ""
		temp_species_df = temp_species_df.loc[species_df["end"] >= row["snp"]]
		if len(temp_species_df) == 0:
			return ""
		temp_species_df = temp_species_df.loc[temp_species_df["protein_id"] != ""]
		if len(temp_species_df) == 0:
			return ""
		else:
			return temp_species_df.iloc[0]["protein_id"]
	else:
		return ""

# identify produced protein
def add_product(row, FP_feature_list):
	if row["species"] in FP_feature_list.keys():
		species_df = FP_feature_list[row["species"]]
		if len(species_df) == 0:
			return ""
		temp_species_df = species_df.loc[species_df["start"] <= row["snp"]]
		if len(temp_species_df) == 0:
			return ""
		temp_species_df = temp_species_df.loc[species_df["end"] >= row["snp"]]
		if len(temp_species_df) == 0:
			return ""
		temp_species_df = temp_species_df.loc[temp_species_df["product"] != ""]
		if len(temp_species_df) == 0:
			return ""
		else:
			return temp_species_df.iloc[0]["product"]
	else:
		return ""


# read in FP_feature
FP_feature_list=defaultdict()
FP_feature_dir="FP_feature"
for f in os.listdir(FP_feature_dir):
	if f.endswith(".gff3"):
		df = pd.read_csv(os.path.join(FP_feature_dir, f), skiprows=2,
			sep='\t', header=None, names=["id", "seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"])
		df = df.drop(["id", "score", "strand", "phase"], axis=1)
		df = df.loc[df["type"] != "region"]
		df["id"] = df["attributes"].map(lambda x: x.split('ID=')[1].split(';')[0] if "ID=" in x else "")
		df["protein_id"] = df["attributes"].map(lambda x: x.split('protein_id=')[1].split(';')[0] if "protein_id=" in x else x.split("RefSeq:")[1].split(';')[0] if "RefSeq:" in x else "")
		df["product"] = df["attributes"].map(lambda x: x.split('product=')[1].split(';')[0] if "product=" in x else "")
		FP_feature_list[f[:-5].replace('.', '_')] = df

# read in SNPs
SNP_df = pd.read_csv("FP_2.csv")

# add columns
SNP_df["gene_id"] = SNP_df.apply(lambda row: add_gene_id(row, FP_feature_list), axis=1)
SNP_df["start"] = SNP_df.apply(lambda row: add_start(row, FP_feature_list), axis=1)
SNP_df["end"] = SNP_df.apply(lambda row: add_end(row, FP_feature_list), axis=1)
SNP_df["protein_id"] = SNP_df.apply(lambda row: add_protein_id(row, FP_feature_list), axis=1)
SNP_df["product"] = SNP_df.apply(lambda row: add_product(row, FP_feature_list), axis=1)

# output mapped gene list
SNP_df.to_csv(sys.argv[3], index=False)

# output mapped gene acc
SNP_acc_df = SNP_df.drop(["Sample", "snp"], axis=1)
SNP_acc_df = SNP_acc_df.groupby(["species", "gene_id", "start", "end", "protein_id", "product"]).sum().reset_index()
SNP_acc_df.to_csv(sys.argv[4], index=False)