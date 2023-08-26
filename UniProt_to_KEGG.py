#!/usr/bin/env python3

name = 'UniProt_to_KEGG.py'
updated = '2023-08-12'
version = '0.5.0'

def UniProt_to_KEGG(parsed_peptide_dir=None,ref=None,upa_seq=None,kegg_id=None,outdir="UniProt_to_KEGG",upa_ecn=None):

	from os import makedirs, listdir
	from os.path import isdir, isfile
	from shutil import which
	from subprocess import run

	## Check for BLASTP binary, kill the process if there isn't one
	blast_type = False
	if which("diamond"):
		blast_type = "DIAMOND"

	elif which("blastp"):
		blast_type = "NCBI"

	if not blast_type:
		print("\n  [E] Could not find a 'blastp' program on the PATH! Please add either NCBI or DIAMOND blastp to your PATH. Exiting...\n")
		exit()


	if not isdir(outdir):
		makedirs(outdir,mode=0o755)
	

	## Run a BLASTP of UniProt sequences against the NCBI reference to link the UniProt accessions to the NCBI accessions
	if not isfile(f"{outdir}/UniProt_to_NCBI.blastp.6"):

		print("\n\tRunning BLASTP search...")

		if blast_type == "DIAMOND":
			
			print("\n\tMaking DIAMOND blastp")
			run(["diamond","makedb",
				"--in",ref,
				"--db",f"{outdir}/diamond_db"],
				capture_output=True)
			
			print("\tRunning DIAMOND blastp")
			run(["diamond","blastp",
				"--query",upa_seq,
				"--db",f"{outdir}/diamond_db",
				"--outfmt","6",
				"qseqid","sseqid",
				"--out",f"{outdir}/UniProt_to_NCBI.blastp.6",
				"--max-target-seqs","1",
				"--max-hsps","1"],
				capture_output=True)

		elif blast_type == "NCBI":

			print("\n\tRunning NCBI blastp")
			run(["blastp",
				"-query",upa_seq,
				"-subject",ref,
				"-outfmt",'6 qseqid sseqid',
				"-out",f"{outdir}/UniProt_to_NCBI.blastp.6",
				"-culling_limit","1",
				"-max_hsps","1"],
				capture_output=True)


	if not isfile(f"{outdir}/UniProt_to_NCBI.blastp.6"):
		print(f"\n  [E] {blast_type} blastp did not execute as expected... Exiting...\n")
		exit()


	## Link EC numbers to UniProt accession numbers
	print("\n\tLoading in Enzyme Classification Numbers...")
	ecns = {}
	ECN = open(upa_ecn,'r')
	for line in ECN:
		line = line.strip()
		upa,ecn = line.split("\t")
		if ecn != "N/A":
			ecns[upa] = ecn
	ECN.close()

	## Getting NCBI accession numbers
	print("\n\tGetting reference accession numbers...")
	accessions = {}
	REF = open(ref,'r')
	for line in REF:
		line = line.strip()
		if line != "":
			if line[0] == ">":
				accessions[line[1:].split(" ")[0].split(".")[0]] = {}
	REF.close()

	## Link NCBI accession numbers to UniProt accessions
	print("\n\tLinking UniProt and NCBI acession numbers...")
	upa_to_ncbi = {}
	ncbi_to_upa = {}
	BLAST = open(f"{outdir}/UniProt_to_NCBI.blastp.6",'r')
	for line in BLAST:
		
		line = line.strip()
		qseqid,sseqid = line.split("\t")[0:2]
		sseqid = sseqid.split(".")[0]

		upa_to_ncbi[qseqid] = sseqid

		if qseqid in ecns.keys():
			accessions[sseqid]['ecn'] = ecns[qseqid]

		if sseqid not in ncbi_to_upa.keys():
			ncbi_to_upa[sseqid] = []

		ncbi_to_upa[sseqid].append(qseqid)

	BLAST.close()

	## Link NCBI accessions to KEGG IDs
	if not isfile(f"{outdir}/{kegg_id}.conv"):
		run(["curl",f"https://rest.kegg.jp/conv/{kegg_id}/ncbi-proteinid","-o",f"{outdir}/{kegg_id}.conv"],capture_output=True)

	kid_to_ncbi = {}
	CONV = open(f"{outdir}/{kegg_id}.conv",'r')
	for line in CONV:
		line = line.strip()
		if line != "":
			ncbi,kid = line.split("\t")
			kid_to_ncbi[kid.replace(f"{kegg_id}:","")] = ncbi.replace(f"ncbi-proteinid:","")
			accessions[ncbi.replace(f"ncbi-proteinid:","")]['kid'] = kid.replace(f"{kegg_id}:","")
	CONV.close()


	## Link NCBI accessions to KEGG Orthologs
	if not isfile(f"{outdir}/{kegg_id}.ko"):
		run(["curl",f"https://rest.kegg.jp/link/ko/{kegg_id}","-o",f"{outdir}/{kegg_id}.ko"],capture_output=True)

	print("\n\tLinking KEGG Ortholog and KEGG IDs")
	KOS = open(f"{outdir}/{kegg_id}.ko",'r')
	for line in KOS:
		
		line = line.strip()
		
		if line != "":

			kid,ko = line.split("\t")
			kid = kid.replace(f"{kegg_id}:","").split(".")[0]
			ko = ko.replace("ko:","")

			if kid in kid_to_ncbi.keys():
				accessions[kid_to_ncbi[kid]]['ko'] = ko

	KOS.close()


	## Link NCBI to KEGG_ID, KEGG_Ortholog_ID, and EC_Number
	print("\n\tLinking NCBI accessions to KEGG Orthologs, KEGG IDs, and Enzyme Classification Numbers...")
	OUT = open(f"{outdir}/conversions.tsv",'w')
	for ref_access in sorted(accessions.keys()):

		OUT.write(f"{ref_access}")

		if ref_access in ncbi_to_upa.keys():
			upas = ";".join(ncbi_to_upa[ref_access])
			OUT.write(f"\t{upas}")
		else:
			OUT.write("\tN/A")


		if 'kid' in accessions[ref_access].keys():
			OUT.write(f"\t{accessions[ref_access]['kid']}")
		else:
			OUT.write("\tN/A")


		if 'ko' in accessions[ref_access].keys():
			OUT.write(f"\t{accessions[ref_access]['ko']}")
		else:
			OUT.write(f"\tN/A")


		if 'ecn' in accessions[ref_access].keys():
			OUT.write(f"\t{accessions[ref_access]['ecn']}")
		else:
			OUT.write(f"\tN/A")


		OUT.write("\n")

	OUT.close()


	print("\n\tConverting Sample accessions to NCBI accessions...")
	sample_ids = [i for i in listdir(parsed_peptide_dir) if isdir(f"{parsed_peptide_dir}/{i}")]
	accession_by_sample = {x:{y:0 for y in sample_ids} for x in accessions.keys()}

	for item in sorted(listdir(parsed_peptide_dir)):

		if isdir(f"{parsed_peptide_dir}/{item}"):
			
			UPA = open(f"{parsed_peptide_dir}/{item}/accessions.list",'r')

			for line in UPA:

				line = line.strip()

				if line in upa_to_ncbi:

					accession_by_sample[upa_to_ncbi[line]][item] = 1

			UPA.close()

	OUT = open(f"{outdir}/accession_by_sample.tsv",'w')
	for accession in sorted(accession_by_sample.keys()):

		OUT.write(f"{accession}")
		
		for item in sorted(accession_by_sample[accession].keys()):

			OUT.write(f"\t{item}:{accession_by_sample[accession][item]}")

		OUT.write("\n")

	OUT.close()

if __name__ == "__main__":


	from argparse import ArgumentParser

	GetOptions = ArgumentParser()
	GetOptions.add_argument("-p","--par_pep_dir",required=True)
	GetOptions.add_argument("-s","--upa_seq",required=True)
	GetOptions.add_argument("-e","--upa_ecn",default=None)
	GetOptions.add_argument("-r","--ref_faa",required=True)
	GetOptions.add_argument("-k","--kegg_id",required=True)
	GetOptions.add_argument("-o","--outdir",default="UniProt_to_KEGG")

	args = GetOptions.parse_known_args()[0]

	peptide_dir = args.par_pep_dir
	upa_seq = args.upa_seq
	ref = args.ref_faa
	upa_ecn = args.upa_ecn
	kegg_id = args.kegg_id
	outdir = args.outdir

	UniProt_to_KEGG(parsed_peptide_dir=peptide_dir,upa_seq=upa_seq,ref=ref,kegg_id=kegg_id,upa_ecn=upa_ecn,outdir=outdir)
