#!/usr/bin/env python3

name = 'make_master_table.py'
updated = '2023-08-31'
version = '0.4.0'

def make_master_table(prot_count_dir=None,ref_faa=None,conversion_file=None,kegg_dir=None,outdir=None):
	
	from os import listdir, makedirs
	from os.path import isdir

	if not isdir(outdir):
		makedirs(outdir,mode=0o755)

	annotations = {}
	table_data = {}
	accession_data = {}
	ANNOT = open(ref_faa,'r')
	for line in ANNOT:
		line = line.strip()
		if line != "":
			if line[0] == ">":
				accession = line[1:].split(" ")[0].split(".")[0]
				annotation = " ".join(line.split(" ")[1:])
				annotations[accession] = annotation
				table_data[accession] = {'annot':annotation}
				accession_data[accession] = []
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


	path_by_ortho = {}
	ORTH = open(f"{kegg_dir}/etc/ortho_pathways.tsv",'r')
	for line in ORTH:
		data = line.strip().split("\t")
		ortho = data[0]
		pathways = []
		for item in data[1].split(";"):
			pathways.append(item.split(":")[0])
		path_by_ortho[ortho] = ";".join(pathways)
	ORTH.close()


	headers = False
	sample_data = {}
	for sample in listdir(prot_count_dir):

		if sample != "temp_dir" and isdir(f"{prot_count_dir}/{sample}"):

			sample_data[sample] = []

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
							if i == 0:
								if int(x) > 0:
									accession_data[accession].append(sample)
									sample_data[sample].append(accession)
							if headers[i] not in table_data[accession].keys():
								table_data[accession][headers[i]] = {}
							table_data[accession][headers[i]][sample] = x

			COUNT.close()

	exclusivity_data = {}
	sample_ids = [x for x in sorted(sample_data.keys())]
	for accession in accession_data.keys():
		exc_string = ""
		for sam in sample_ids:
			if sam in accession_data[accession]:
				exc_string = '1' + exc_string
			else:
				exc_string = '0' + exc_string
		if exc_string not in exclusivity_data.keys():
			exclusivity_data[exc_string] = []
		exclusivity_data[exc_string].append(accession)

	exc_ids = [x for x in sorted(exclusivity_data.keys())]
	exc_dataset	= [exclusivity_data[x] for x in exc_ids]
	max_iter = max([len(x) for x in exc_dataset])
	
	OUT = open(f"{outdir}/ExclusivityGroups.tsv",'w')
	title1 = "\t".join(exc_ids)
	OUT.write(f"## {title1}\n")
	title2 = [":".join([sample_ids[len(sample_ids)-index-1] for index,char in enumerate(exc) if char == '1']) for exc in exc_ids]
	title2[0] = 'MissingInAllSamples'
	title2 = "\t".join(title2)
	OUT.write(f"# {title2}\n")
	for index in range(max_iter):
		line = []
		for dataset in exc_dataset:
			if index < len(dataset):
				line.append(dataset[index])
			else:
				line.append("")
		line = "\t".join(line)
		OUT.write(f"{line}\n")
	OUT.close()


	def write_master(outpath,accessions_list,sample,samples=None):

		OUT = open(outpath,'w')
		
		if sample == 'all':
			header_1 = ("\t"*len(sample_ids)).join(headers)
			header_2 = "\t".join(("\t".join(sample_ids)) for x in range(len(headers)))
			OUT.write(f"## Accession\tProvided_UniProt_Accessions\tEnzyme_Classification_Number\tKEGG_ID\tKEGG_Ortholog\tPathways\tAnnotation\tSequence_Length\t{header_1}\n")
			OUT.write(f"#\t\t\t\t\t\t\t\t{header_2}\n")
		elif sample == "exc":
			header_1 = ("\t"*len(samples)).join(headers)
			header_2 = "\t".join(("\t".join(samples)) for x in range(len(headers)))
			OUT.write(f"## Accession\tProvided_UniProt_Accessions\tEnzyme_Classification_Number\tKEGG_ID\tKEGG_Ortholog\tPathways\tAnnotation\tSequence_Length\t{header_1}\n")
			OUT.write(f"#\t\t\t\t\t\t\t\t{header_2}\n")
		else:
			header_1 = ("\t").join(headers)
			OUT.write(f"## Accession\tProvided_UniProt_Accessions\tEnzyme_Classification_Number\tKEGG_ID\tKEGG_Ortholog\tPathways\tAnnotation\tSequence_Length\t{header_1}\n")

		for accession in accessions_list:

			OUT.write(accession)

			upas = table_data[accession]['upas']
			slen = table_data[accession]['len']
			annot = table_data[accession]['annot']
			kid = table_data[accession]['kid']
			ko = table_data[accession]['ko']
			ecn = table_data[accession]['ecn']
			
			pathways = ""
			if ko in path_by_ortho.keys():
				pathways = path_by_ortho[ko]
			else:
				pathways = "N/A"

			OUT.write(f"\t{upas}\t{ecn}\t{kid}\t{ko}\t{pathways}\t{annot}\t{slen}")

			for head in headers:

				if sample == 'all':

					for sam in sample_ids:

						OUT.write(f"\t{table_data[accession][head][sam]}")

				elif sample == 'exc':

					for sam in samples:

						OUT.write(f"\t{table_data[accession][head][sam]}")
				
				else:

					OUT.write(f"\t{table_data[accession][head][sample]}")

			OUT.write("\n")

		OUT.close()

	write_master(f"{outdir}/master_table.tsv",sorted(table_data.keys()),'all')

	temp_dir = f"{outdir}/sample_tables"
	if not isdir(temp_dir):
		makedirs(temp_dir,mode=0o755)
	for sam in sample_ids:
		write_master(f"{temp_dir}/{sam}_master_table.tsv",sample_data[sam],sam)
	for exc in exc_ids:
		write_master(f"{temp_dir}/{exc}_master_table.tsv",exclusivity_data[exc],'exc',samples=[sample_ids[i] for i,x in enumerate(reversed(exc)) if x != '0'])
	
	
	OUT = open(f"{outdir}/accession_per_sample.tsv",'w')
	max_iter = max([len(sample_data[x]) for x in sample_data.keys()])
	headers = '# ' + "\t\t".join([x for x in sorted(sample_data.keys())])
	OUT.write(f"{headers}\n")
	for i in range(max_iter):
		line = []
		for sam in sorted(sample_data.keys()):
			if i < len(sample_data[sam]):
				line.append(sample_data[sam][i])
				line.append(annotations[sample_data[sam][i]])
			else:
				line.append("")
				line.append("")
		line = "\t".join(line)
		OUT.write(f"{line}\n")
	OUT.close()

	
	OUT = open(f"{outdir}/simple_exclusivity_table.tsv",'w')
	title1 = "\t\t\t".join(exc_ids)
	OUT.write(f"### {title1}\n")
	title2 = [":".join([sample_ids[len(sample_ids)-index-1] for index,char in enumerate(exc) if char == '1']) for exc in exc_ids]
	title2[0] = 'MissingInAllSamples'
	title2 = "\t\t\t".join(title2)
	OUT.write(f"## {title2}\n")
	title3 = "\t".join([f"NCBI Accession\tAnnotation\tPathways" for x in range(len(exc_ids))])
	OUT.write(f"# {title3}\n")
	for index in range(max_iter):
		line = []
		for dataset in exc_dataset:
			if index < len(dataset):
				accession = dataset[index]
				line.append(accession)
				annot = table_data[accession]['annot']
				line.append(annot)
				ko = table_data[accession]['ko']
				if ko in path_by_ortho.keys():
					line.append(path_by_ortho[ko])
				else:
					line.append("N/A")
			else:
				line.append("")
				line.append("")
				line.append("")
		line = "\t".join(line)
		OUT.write(f"{line}\n")
	OUT.close()

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()

	GetOptions.add_argument("-p","--prot_count_dir",required=True)
	GetOptions.add_argument("-r","--ref_faa",required=True)
	GetOptions.add_argument("-c","--conv_file",required=True)
	GetOptions.add_argument("-k","--kegg_dir",required=True)
	GetOptions.add_argument("-o","--outdir",default="SummarizedData")

	args = GetOptions.parse_known_args()[0]
	prot_count_dir = args.prot_count_dir
	conv_file = args.conv_file
	ref_faa = args.ref_faa
	kegg_dir = args.kegg_dir
	outdir = args.outdir


	make_master_table(prot_count_dir=prot_count_dir,ref_faa=ref_faa,conversion_file=conv_file,kegg_dir=kegg_dir,outdir=outdir)