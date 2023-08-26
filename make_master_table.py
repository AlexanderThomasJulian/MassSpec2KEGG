#!/usr/bin/env python3

name = 'make_master_table.py'
updated = '2023-08-12'
version = '0.2.0'

def make_master_table(prot_count_dir=None,ref_faa=None,conversion_file=None):
	
	from os import listdir
	from os.path import isdir

	annotations = {}
	table_data = {}
	ANNOT = open(ref_faa,'r')
	for line in ANNOT:
		line = line.strip()
		if line != "":
			if line[0] == ">":
				accession = line[1:].split(" ")[0].split(".")[0]
				annotation = " ".join(line.split(" ")[1:])
				annotations[accession] = annotation
				table_data[accession] = {'annot':annotation}
	ANNOT.close()


	CONV = open(conversion_file,'r')
	for line in CONV:
		line = line.strip()
		accession,upas,kid,ko,ecn = line.split("\t")
		table_data[accession]['upas'] = upas
		table_data[accession]['kid'] = kid
		table_data[accession]['ko'] = ko
		table_data[accession]['ecn'] = ecn
	CONV.close()

	headers = False
	samples = {}
	for sample in listdir(prot_count_dir):

		if sample != "temp_dir" and isdir(f"{prot_count_dir}/{sample}"):

			samples[sample]= True

			COUNT = open(f"{prot_count_dir}/{sample}/count.tsv",'r')

			for line in COUNT:
				line = line.strip()
				if line != "":
					if line[0] == "#":
						if not headers:
							headers = line.split("\t")[2:]
					else:
						data = line.split("\t")
						accession = data[0]
						seq_len = data[1]

						table_data[accession]['len'] = seq_len
						for i,x in enumerate(data[2:]):
							if headers[i] not in table_data[accession].keys():
								table_data[accession][headers[i]] = {}
							table_data[accession][headers[i]][sample] = x

			COUNT.close()


	OUT = open(f"{prot_count_dir}/master_table.tsv",'w')

	
	header_1 = ("\t"*len(samples.keys())).join(headers)
	header_2 = "\t".join(("\t".join(sorted(samples.keys()))) for x in range(len(headers)))
	OUT.write(f"## Accession\tProvided_UniProt_Accessions\tKEGG_ID\tKEGG_Ortholog\tEnzyme_Classification_Number\tAnnotation\tSequence_Length\t{header_1}\n")
	OUT.write(f"##\t\t\t\t\t\t\t{header_2}\n")

	for accession in table_data.keys():

		OUT.write(accession)

		upas = table_data[accession]['upas']
		slen = table_data[accession]['len']
		annot = table_data[accession]['annot']
		kid = table_data[accession]['kid']
		ko = table_data[accession]['ko']
		ecn = table_data[accession]['ecn']

		OUT.write(f"\t{upas}\t{kid}\t{ko}\t{ecn}\t{annot}\t{slen}")

		for head in headers:

			for sample in sorted(table_data[accession][head].keys()):

				OUT.write(f"\t{table_data[accession][head][sample]}")

		OUT.write("\n")

	OUT.close()


if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-p","--prot_count_dir",required=True)
	GetOptions.add_argument("-r","--ref_faa",required=True)
	GetOptions.add_argument("-c","--conv_file",required=True)

	args = GetOptions.parse_known_args()[0]
	prot_count_dir = args.prot_count_dir
	conv_file = args.conv_file
	ref_faa = args.ref_faa


	make_master_table(prot_count_dir=prot_count_dir,ref_faa=ref_faa,conversion_file=conv_file)