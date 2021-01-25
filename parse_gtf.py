#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import re
from argparse import ArgumentParser
import sys
import glob
import os

args = ArgumentParser('./parse_gtf.py', description="""This program has been designed to obtain
information about transcripts from a GTF file. A .csv file containing information about each transcript
will be generated. In addition, if you would like a separate file with the information only about the
longest transcript per gene or the transcript with the longest coding sequence, you can specify to have
those files created with the -t or -c options.
Example usage: python parse_gtf.py --gtf_file Homo_sapiens.GRCh38.96.gtf
--output_prefix HG38.96 -c""")

args.add_argument(
	'-g',
	'--gtf_file',
	help="""This is the GTF file from which you would like to parse information.
	If you do not provide a file, the program will attempt to find a GTF file in your present working directory.""",
	default=None,
)

args.add_argument(
	'-o',
	'--output_prefix',
	help="This is the name of the output file that you would like created. If you do not specify a prefix, the prefix will be taken from the GTF file that you provided.",
	default = None
)

args.add_argument(
	'-t',
	'--max_transcript_length',
	help="Specify this option if you would also like to make a file with only the transcripts containing the longest transcript length.",
	action= 'store_true'
)

args.add_argument(
	'-c',
	'--max_cds',
	help="Specify this option if you would also like to make a file with only the transcripts containing the longest coding sequence (CDS).",
	action= 'store_true'
)

args = args.parse_args()
input_file = args.gtf_file
if not input_file:
	gtf_files = glob.glob('*.gtf')
	if len(gtf_files) == 0:
		print()
		print("Welcome to parse_gtf.py.")
		print("This program has been written to parse out information about each transcript from an input GTF file.")
		print("No GTF files are present in the present working directory.")
		print("Please use the --gtf_file option to specify your path to the GTF file you wish to parse.")
		print("Example usage: python parse_gtf.py --gtf_file Homo_sapiens.GRCh38.96.gtf --output_prefix HG38.96 -c")
		print()
		sys.exit(1) # Exit with a status of 1. They are probably trying to see what the program does.
	elif len(gtf_files) > 1:
		print()
		print("Welcome to parse_gtf.py.")
		print("This program has been written to parse out information about each transcript from an input GTF file.")
		print("You have multiple GTF files present in your working directory.")
		print("Please specify a GTF file with the --gtf_file option.")
		print("The GTF files in your current working directory are:")
		print(gtf_files)
		print("Example usage: python parse_gtf.py --gtf_file Homo_sapiens.GRCh38.96.gtf --output_prefix HG38.96 -c")
		print()
		sys.exit(1)
	elif len (gtf_files) == 1: # This should be all others, but just to be safe I am using ==1
		input_file = gtf_files[0]

prefix = args.output_prefix
if not prefix:
	prefix = os.path.splitext(input_file)[0]
transcript_output_file = prefix+"_Transcript_information.csv"
max_transcript_length = args.max_transcript_length
max_cds = args.max_cds

count = 0
for i, line in enumerate(open(input_file)):
	if line.startswith('#!'):
		count += 1
	else:
		break

gtf = pd.read_csv(input_file, skiprows = count, index_col = None, header = None, sep = '\t',
					dtype = {0: str, 1: str, 2: str, 3: int, 4: int, 5: str, 6: str,7: str,8: str})
gtf.columns = ["Chr", "Source", "Feature", "Start", "Stop",  "Score", "Strand", "Frame", "Info"]

transcripts = gtf[gtf['Feature'] == 'transcript']
def parse_transcript_info(df):
	info = df['Info']
	start = df['Start']
	stop = df['Stop']
	length = (stop-start)+1
	gene_id = None
	gene_version = None
	transcript_id = None
	transcript_version = None
	gene_name = None
	# gene_source # This one is already present
	gene_biotype = None
	transcript_name = None
	# transcript_source # I am just going to use the gene source
	tag = None
	transcript_support_level = None
	transcript_biotype = None
	### Now I need to split the info field into a dictionary and put the into the dataframe
	info_dict = dict(x.split(' "') for x in info.split("; "))
	### Now update the variables created above if they are present in the info field
	variable_info = {
		'gene_id': gene_id,
		'gene_version': gene_version,
		'transcript_id': transcript_id,
		'transcript_version': transcript_version,
		'gene_name': gene_name,
		'gene_biotype': gene_biotype,
		'transcript_name': transcript_name,
		'tag': tag,
		'transcript_support_level': transcript_support_level,
		'transcript_biotype': transcript_biotype
	}
	list_of_variables = [gene_id, gene_version, transcript_id, transcript_version, gene_name, gene_biotype, transcript_name, tag, transcript_support_level, transcript_biotype]
	list_of_strings = ["gene_id", "gene_version", "transcript_id", "transcript_version", "gene_name", "gene_biotype", "transcript_name", "tag", "transcript_support_level", "transcript_biotype"]
	for string, variable in variable_info.items():
		if string in info_dict:
			variable = info_dict[string]
			variable = variable.replace('"', '').replace(';', '')
		df[string] = variable
	df['Transcript Length'] = length
	return df


transcripts = transcripts.apply(parse_transcript_info, axis = 1)

exons = gtf[gtf['Feature'] == 'exon']
def parse_exon_info(df):
	info = df['Info']
	start = df['Start']
	stop = df['Stop']
	length = (stop - start)+1
	transcript_id = None
	info_dict = dict(x.split(' "') for x in info.split("; "))
	if "transcript_id" in info_dict:
		transcript_id = info_dict['transcript_id'].replace('"', '').replace(';', '')
	df['CDS Length'] = length
	df['transcript_id'] = transcript_id
	return df


exons = exons.apply(parse_exon_info, axis = 1)
sums = exons.groupby('transcript_id').sum()
exonic_lengths = sums[['CDS Length']]
transcripts = transcripts.merge(exonic_lengths, on = 'transcript_id')

def get_intronic_length(df):
	total_len = df['Transcript Length']
	exon_len = df['CDS Length']
	intronic_len = total_len-exon_len
	df['Intronic Length'] = intronic_len
	return df


transcripts = transcripts.apply(get_intronic_length, axis = 1)
transcripts = transcripts[['transcript_id', 'gene_id', 'Chr', 'Start', 'Stop',
							'Strand', 'Transcript Length','CDS Length', 'Intronic Length',
							'gene_name', 'gene_biotype', 'gene_version','transcript_name',
							'transcript_biotype', 'transcript_version', 'transcript_support_level','tag', 'Source']]

transcripts.to_csv(transcript_output_file)

## Now to make the files with the longest transcript or the longest CDS length

def find_longest(df, column_name, lengths_dict):
	gene_id = df['gene_id']
	length = df[column_name]
	longest = lengths_dict[gene_id]
	keep = False # set the default value to not keep the row
	if length == longest:
		keep = True
	df['Longest'] = keep
	return df


if max_transcript_length:
	longest_transcript_output = prefix+"_Longest_Transcript.csv"
	lengths_df = pd.DataFrame(transcripts.groupby(['gene_id'])['Transcript Length'].max())
	transcript_lengths_dict = dict(zip(lengths_df.index.values.tolist(), lengths_df['Transcript Length'].values.tolist()))
	longest_transcript_df = transcripts.apply(find_longest, column_name = "Transcript Length", lengths_dict = transcript_lengths_dict, axis = 1)
	longest_transcript_df = longest_transcript_df[longest_transcript_df['Longest'] == True]
	longest_transcript_df = longest_transcript_df[longest_transcript_df.columns[:-1]]
	### Some of them have multiple transcripts with the same length. I am going to just take the first entry for those ones.
	gene_ids = longest_transcript_df["gene_id"]
	bool_duplicates = dict(gene_ids.duplicated("first"))
	duplicate_results = {}
	for entry in bool_duplicates:
		if bool_duplicates[entry]: # If it comes back as true
			duplicate_results[entry] = bool_duplicates[entry]
	dup_indexes = list(duplicate_results.keys()) # This is the indexes of the ones that returned "True"
	updated_df = longest_transcript_df.drop(dup_indexes)
	updated_df.to_csv(longest_transcript_output)

if max_cds:
	longest_coding_output = prefix+"_Longest_CDS.csv"
	lengths_df = pd.DataFrame(transcripts.groupby(['gene_id'])['CDS Length'].max())
	exonic_lengths_dict = dict(zip(lengths_df.index.values.tolist(), lengths_df['CDS Length'].values.tolist()))
	longest_exonic_df = transcripts.apply(find_longest, column_name = "CDS Length", lengths_dict = exonic_lengths_dict, axis = 1)
	longest_exonic_df = longest_exonic_df[longest_exonic_df['Longest'] == True]
	longest_exonic_df = longest_exonic_df[longest_exonic_df.columns[:-1]]
	### Some of them have multiple transcripts with the same CDS length. I am going to just take the first entry for those ones.
	gene_ids = longest_exonic_df["gene_id"]
	bool_duplicates = dict(gene_ids.duplicated("first"))
	duplicate_results = {}
	for entry in bool_duplicates:
		if bool_duplicates[entry]: # If it comes back as true
			duplicate_results[entry] = bool_duplicates[entry]
	dup_indexes = list(duplicate_results.keys()) # This is the indexes of the ones that returned "True"
	updated_df = longest_exonic_df.drop(dup_indexes)
	updated_df.to_csv(longest_coding_output)
