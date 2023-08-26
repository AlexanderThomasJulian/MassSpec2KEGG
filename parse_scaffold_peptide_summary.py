#!/usr/bin/env python3

name = 'parse_scaffolds_peptide_summary.py'
updated = '2023-08-02'
version = '0.1.0'

def parse_scaffold_peptide_summary(peptide_report=False,outdir="PARSED_REPORT"):

	from os.path import isdir,isfile
	from os import makedirs

	## Verify that file is provided
	if not peptide_report or not isfile(peptide_report):
		print(f"\n  [E] peptide_report ({peptide_report}) is not accessible. Be sure the file provided exists. Exiting...\n")
		exit()


	grab = False
	REPORT = open(report,'r')
	data = {}
	accessions = {}
	sample = ""
	for count,line in enumerate(REPORT):
		line = line.strip()

		if grab:
			if line != "" and len(line.split("\t")) > 1:
				
				temp = line.split("\t")
				sample = temp[1]
				accession = temp[5].split(",")[0]
				id_prob = float(temp[9].replace("%",""))
				seq = temp[15]

				if id_prob > 80:

					accessions[accession] = True
					
					if sample not in data.keys():
						data[sample] = {}
					
					if accession not in data[sample].keys():
						data[sample][accession] = []
					
					data[sample][accession].append(seq)
			
		if grab == False and len(line.split("\t")) == 37:
			grab = True

	REPORT.close()


	MASTER_ACCESSION = open(f"{outdir}/accessions.list",'w')
	for accession in sorted(accessions.keys()):
		MASTER_ACCESSION.write(f"{accession}\n")
	MASTER_ACCESSION.close()


	if not isdir(outdir):
		makedirs(outdir,mode=0o755)

	for sample in data.keys():
		
		sample_name = sample.replace(" ","_")
		temp_dir = f"{outdir}/{sample_name}"
		
		if not isdir(temp_dir):
			makedirs(temp_dir,mode=0o755)
		
		if not isdir(f"{temp_dir}/sequences"):
			makedirs(f"{temp_dir}/sequences",mode=0o755)
		
		ACCESSIONS = open(f"{temp_dir}/accessions.list",'w')
		
		for accession in sorted(data[sample].keys()):

			ACCESSIONS.write(f"{accession}\n")

			SEQ = open(f"{temp_dir}/sequences/{accession}.faa",'w')

			for count,seq in enumerate(sorted(data[sample][accession])):

				SEQ.write(f">{accession}:{count+1}\n{seq}\n")

			SEQ.close()

		ACCESSIONS.close()

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()
	GetOptions.add_argument("-p","--peptide_report",required=True)
	GetOptions.add_argument("-o","--outdir",default="PARSED_REPORT")

	args = GetOptions.parse_known_args()[0]
	report = args.peptide_report
	outdir = args.outdir

	parse_scaffold_peptide_summary(report,outdir)