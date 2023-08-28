#!/usr/bin/env python3

name = 'protein_count.py'
version = '0.1.1'
updatd = '2023-08-28'

from argparse import ArgumentParser
from os.path import isdir, isfile
from os import listdir, makedirs
from subprocess import run
from progress_bar import ProgressBar
from time import time

GetOptions = ArgumentParser()
GetOptions.add_argument("-b","--utn_blast",required=True)
GetOptions.add_argument("-n","--ncbi_faa",required=True)
GetOptions.add_argument("-p","--parsed_dir",required=True)
GetOptions.add_argument("-o","--outdir",default="Protein_Count")

args = GetOptions.parse_known_args()[0]
blast = args.utn_blast
ref = args.ncbi_faa
parsed_dir = args.parsed_dir
outdir = args.outdir

if not isdir(outdir):
	makedirs(outdir,0o755)

if not isdir(f"{outdir}/temp_dir"):
	makedirs(f"{outdir}/temp_dir",0o755)

if not isdir(f"{outdir}/temp_dir/split_ref"):
	makedirs(f"{outdir}/temp_dir/split_ref",0o755)

accessions = {}
annotations = {}
locus = ""

## Load in the NCBI reference FASTA, and store the sequences
REF = open(ref,'r')
print(f"\tLoading reference...")
for line in REF:
	line = line.strip()
	if line != "":
		if line[0] == ">":
			locus = line[1:].split(" ")[0].split(".")[0]
			accessions[locus] = {
				'seq':"",
				'len':0,
				'matches':[]
			}
		else:
			accessions[locus]['seq'] += line
			accessions[locus]['len'] += len(line)
REF.close()

## Store each reference sequence in its own fasta file. Needed to perform piece-wise BLAST search. If skipped, may return
## spurious hits that do not match the SCAFFOLD inference
print(f"\tSplitting reference...")
for count,locus in enumerate(accessions.keys()):
	if not isfile(f"{outdir}/temp_dir/split_ref/{locus}.faa"):
		OUT = open(f"{outdir}/temp_dir/split_ref/{locus}.faa",'w')
		OUT.write(f">{locus}\n{accessions[locus]}\n")
		OUT.close()


## Read in the links between UNIPROT accessions and NCBI accessions
BLAST = open(blast,'r')
print(f"\tReading BLAST...")
total_matches = 0
for line in BLAST:
	
	line = line.strip()
	
	uniprot,ncbi = line.split("\t")
	ncbi = ncbi.split(".")[0]
	
	accessions[ncbi]['matches'].append(uniprot)
	total_matches += 1

BLAST.close()


### NEEDS TO BE OPTIMIZED! Multithreading, maybe?
## Perform piece-wise BLAST between all sequences for a given UNIPROT accession and the NCBI accession

### NEEDS TO BE DISCUSSED! SINGLE AMINO ACID MAY NOT BE STRONG ENOUGH TO SUGGEST COUNT!
## Count each hit at each Amino Acid sequence, highest count of a single Amino Acid is is marked as protein count.

for item in listdir(parsed_dir):

	print(f"\tDetermining Protein Count: {item}")

	# pg = ProgressBar(total_matches)

	if isdir(f"{parsed_dir}/{item}"):

		if not isdir(f"{outdir}/{item}"):
			makedirs(f"{outdir}/{item}",mode=0o755)

		COUNT = open(f"{outdir}/{item}/count.tsv",'w')

		if not isdir(f"{outdir}/temp_dir/blast/{item}"):
			makedirs(f"{outdir}/temp_dir/blast/{item}",mode=0o755)

		COUNT.write(f"## Accession\tAccessionSequenceLength\tDBEntriesCovered\tCumulativePartialPeptideSequences\tCummulativePeptides\tAverageSequenceCoverage\tCalculatedProteinCount\n")

		## Iterate over NCBI accessions
		for accession in sorted(accessions.keys()):

			sequence = [0 for x in range(accessions[accession]['len'])]
			accessions_count = 0
			peptide_count = 0
			partial_seqs = 0

			## Iterate over UNIPROT accession linked to current NCBI accession

			for match in accessions[accession]['matches']:

				temp_file = f"{parsed_dir}/{item}/sequences/{match}.faa"

				if isfile(temp_file):
					
					accessions_count += 1

					TEMP = open(temp_file,'r')
					for line in TEMP:
						line = line.strip()
						if line != "":
							if line[0] == ">":
								partial_seqs += 1
							else:
								peptide_count += len(line)
					TEMP.close()


					if not isfile(f"{outdir}/temp_dir/blast/{item}/{match}_{match}.blastp.6"):

						run(["blastp",
								"-query",f"{temp_file}",
								"-subject",f"{outdir}/temp_dir/split_ref/{accession}.faa",
								"-outfmt","6 qseqid sseqid sstart send slen",
								"-out",f"{outdir}/temp_dir/blast/{item}/{match}.blastp.6",
								"-culling_limit","1",
								"-max_hsps",'1'
						])

					BLAST = open(f"{outdir}/temp_dir/blast/{item}/{match}.blastp.6",'r')
					for line in BLAST:
						line = line.strip()
						sseqid,qseqid,start,end,slen = line.split("\t")
						# print(item,accession,match,line,len(sequence))
						for i in range(int(start)-1,int(end)):
							sequence[i] += 1

				# pg.update_progress_bar()

			# COUNT.write(f"## Accession\tAccessionSequenceLength\tDBEntriesCovered\tCumulativePartialPeptideSequences\tCummulativePeptides\tAverageSequenceCoverage\tCalculatedProteinCount\n")

			protein_count = 0
			coverage = 0

			if accessions_count != 0:
				if max(sequence) == 0:
					protein_count = 'Could Not Be Accurately Determined'
					coverage = 'Could Not Be Accurately Determined'
				else:
					protein_count = max(sequence)
					coverage = f"{(sum(sequence)/len(sequence)):.2f}"
			
			# Accession
			COUNT.write(f"{accession}")
			# AccessionSequenceLength
			COUNT.write(f"\t{accessions[accession]['len']}")
			# NumDBEntries
			COUNT.write(f"\t{accessions_count}")
			# CummulativePeptidesSequences
			COUNT.write(f"\t{partial_seqs}")
			# CummulativePeptides
			COUNT.write(f"\t{peptide_count}")
			# AverageSequenceCoverage
			COUNT.write(f"\t{coverage}")
			# CalculatedProteinCount
			COUNT.write(f"\t{protein_count}\n")

		COUNT.close()
	
	# print()