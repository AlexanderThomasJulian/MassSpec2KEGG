#!/usr/bin/env python3

name = 'grab_UniProt_data.py'
updated = '2023-08-02'
vresion = '0.3.0'

def grab_UniProt_data(accessions=False,outdir="UniProtData"):

	from subprocess import run
	from textwrap import wrap
	from json import loads
	from progress_bar import ProgressBar
	from os.path import isdir,isfile
	from os import makedirs


	if not accessions or not isfile(accessions):
		print(f"\n  [E] accessions ({accessions}) is not accessible. Be sure the file provided exists. Exiting...\n")
		exit()

	if not isdir(outdir):
		makedirs(outdir,mode=0o755)


	accessions = {}
	ACCESS = open(access,'r')
	for line in ACCESS:
		line = line.strip()
		accession = line.split("\t")[0]
		if " " not in accession:
			if accession not in accessions.keys():
				accessions[accession] = True
			
	ACCESS.close()


	previous = {}
	if isfile(f"{outdir}/metadata.log"):
		META = open(f"{outdir}/metadata.log",'r')
		for line in META:
			line = line.strip()
			accession,ec,seq = line.split("\t")
			previous[accession] = {'ecNumber':ec,'seq':seq}
		META.close()


	META = open(f"{outdir}/metadata.log",'w')
	ECN = open(f"{outdir}/UPA_ECN.tsv",'w')
	SEQ = open(f"{outdir}/UPA_SEQ.faa",'w')

	PG = ProgressBar(len(accessions.keys()))

	for count,accession in enumerate(accessions.keys()):

		sequence = "N/A"
		ecNumber = "N/A"

		if accession not in previous.keys():

			item = run(["curl",f"https://rest.uniprot.org/uniprotkb/search?query={accession}&fields=ec,sequence"],capture_output=True,text=True).stdout

			data = loads(item)

			if len(data['results']) > 0:
				
				data = data['results'][0]

				if 'sequence' in data.keys():
				
					sequence = data['sequence']['value']
					seq = "\n".join(wrap(sequence,60))
					SEQ.write(f">{accession}\n{seq}\n")

				else:

					item2 = run(["curl",f"https://rest.uniprot.org/uniparc/search?query={accession}&fields=sequence"],capture_output=True,text=True).stdout

					data2 = loads(item2)['results'][0]

					if 'sequence' in data2.keys():

						sequence = data2['sequence']['value']
						seq = "\n".join(wrap(sequence,60))
						SEQ.write(f">{accession}\n{seq}\n")


				if 'proteinDescription' in data.keys():
				
					data = data['proteinDescription']
				
					if 'recommendedName' in data.keys():
				
						data = data['recommendedName']

						if 'ecNumbers' in data.keys():
							ecNumber = data['ecNumbers'][0]['value']
							ECN.write(f"{accession}\t{ecNumber}\n")

		else:

			sequence = previous[accession]['seq']
			ecNumber = previous[accession]['ecNumber']
			seq = "\n".join(wrap(sequence,60))

			SEQ.write(f">{accession}\n{seq}\n")
			ECN.write(f"{accession}\t{ecNumber}\n")

		META.write(f"{accession}\t{ecNumber}\t{sequence}\n")

		PG.update_progress_bar()

	print("\n")

	META.close()
	ECN.close()
	SEQ.close()

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()
	GetOptions.add_argument("-a","--accessions",required=True)
	GetOptions.add_argument("-o","--outdir",default="UniProtData")

	args = GetOptions.parse_known_args()[0]
	access = args.accessions
	outdir = args.outdir

	grab_UniProt_data(accessions=access,outdir=outdir)