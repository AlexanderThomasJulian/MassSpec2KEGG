#!/usr/bin/env python3

if __name__ == "__main__":

	from argparse import ArgumentParser
	from os.path import isdir, isfile, basename, abspath,split
	from os import makedirs, listdir
	from sys import argv
	from parse_scaffold_peptide_summary import parse_scaffold_peptide_summary
	from grab_UniProt_data import grab_UniProt_data
	from UniProt_to_KEGG import UniProt_to_KEGG
	from grab_KEGG_maps import grab_KEGG_maps
	from kgml_parser import kgml_parser
	from kgml_drawer import kgml_drawer


	path = split(abspath(argv[0]))[0]

	GetOptions = ArgumentParser()
	GetOptions.add_argument("-p","--peptide_report",required=True)
	GetOptions.add_argument("-r","--ref_faa",required=True)
	GetOptions.add_argument("-k","--kegg_id",required=True)
	GetOptions.add_argument("-o","--outdir",default="UP_to_KP")


	args = GetOptions.parse_known_args()[0]
	peptide_report = args.peptide_report
	ref = args.ref_faa
	kegg_id = args.kegg_id
	outdir = args.outdir

	if not isdir(outdir):
		makedirs(outdir,mode=0o755)

	PEPTIDE_REPORT_DIR = f"{outdir}/Parsed_Peptide_Report"
	UNIPROT_DATA_DIR = f"{outdir}/UniProt_Data"
	UNIPROT_TO_KEGG_DIR = f"{outdir}/UniProt_to_KEGG"
	KEGG_FILES_DIR = f"{outdir}/KEGG_MAP_FILES"

	## Parse the Scaffolds Peptide_Summary
	parse_scaffold_peptide_summary(
		peptide_report=peptide_report,
		outdir=PEPTIDE_REPORT_DIR
	)

	## Scrap the data from UniProt
	grab_UniProt_data(
		accessions=f"{PEPTIDE_REPORT_DIR}/accessions.list",
		outdir=UNIPROT_DATA_DIR
	)

	## Convert UniProt accessions to NCBI accessions and KEGG IDS
	UniProt_to_KEGG(
		UPA_SEQ=f"{UNIPROT_DATA_DIR}/UPA_SEQ.faa",
		ref_faa=ref,
		KEGG_ID=kegg_id,
		outdir=UNIPROT_TO_KEGG_DIR,
		UPA_ECN=f"{UNIPROT_DATA_DIR}/UPA_ECN.tsv"
	)

	## Download KEGG pathway data
	grab_KEGG_maps(
		outdir=KEGG_FILES_DIR
	)

	## Get KEGG pathway location data ##
	kgml_parser(
		kgml_dir=f"{KEGG_FILES_DIR}/kgmls",
		path_list=f"{KEGG_FILES_DIR}/etc/pathways.list",
		outdir=f"{KEGG_FILES_DIR}/locs"
	)

	for item in listdir(PEPTIDE_REPORT_DIR):

		if isdir(item):

			kgml_drawer(
				conversions_file=f"{UNIPROT_TO_KEGG_DIR}/conversions.tsv",
				KEGG_files_dir=KEGG_FILES_DIR,
				KEGG_pathway_locs_dir=f"{KEGG_FILES_DIR}/locs"
			)