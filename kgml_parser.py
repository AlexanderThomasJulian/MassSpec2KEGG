#!/usr/bin/env python3

name = 'kgml_parser.py'
updated = '2023-08-31'
version = '0.3.1'

def kgml_parser(kegg_dir=None):

	from os.path import isdir, isfile, basename
	from os import makedirs, listdir, stat
	import xml.etree.ElementTree as ET

	outdir = f"{kegg_dir}/loc"

	if not isdir(outdir):
		makedirs(outdir,mode=0o755)

	if not isdir(f"{outdir}/ec"):
		makedirs(f"{outdir}/ec",mode=0o755)

	if not isdir(f"{outdir}/ko"):
		makedirs(f"{outdir}/ko",mode=0o755)


	paths = {}
	PL = open(f"{kegg_dir}/etc/pathways.list",'r')
	for line in PL:
		paths[line.strip()] = True
	PL.close()
	
	## Pathway name to pathway key
	pathways = {}

	## Pathway name per kegg ortholog
	path_by_ortho = {}

	## Number of orthologs per kegg pathway
	ortho_by_path = {}

	## Number of ortholog groups per kegg pathway
	ortho_groups_by_path = {}

	for path in sorted(paths.keys()):

		if isfile(f"{kegg_dir}/kgml/ec/ec{path}.kgml"):

			LOCS = open(f"{outdir}/ec/ec{path}.locs",'w')
			if stat(f"{kegg_dir}/kgml/ec/ec{path}.kgml").st_size != 0:
			
				tree = ET.parse(f"{kegg_dir}/kgml/ec/ec{path}.kgml")

				for entry in tree.getroot():
					if entry.tag == 'entry':
						if entry.attrib['type'] == 'enzyme':
							enzyme = entry.attrib['name'].replace("ec:","")
							for graphic in entry:
								x = 0
								y = 0
								w = 0
								h = 0
								if graphic.attrib['type'] == 'rectangle':
									x = graphic.attrib['x']
									y = graphic.attrib['y']
									w = graphic.attrib['width']
									h = graphic.attrib['height']

								LOCS.write(f"{enzyme}")
								LOCS.write(f"\tx:{x}")
								LOCS.write(f"\ty:{y}")
								LOCS.write(f"\tw:{w}")
								LOCS.write(f"\th:{h}")
								LOCS.write(f"\n")
			LOCS.close()

		
		if isfile(f"{kegg_dir}/kgml/ko/ko{path}.kgml"):

			LOCS = open(f"{outdir}/ko/ko{path}.locs",'w')
			if stat(f"{kegg_dir}/kgml/ko/ko{path}.kgml").st_size != 0:
				tree = ET.parse(f"{kegg_dir}/kgml/ko/ko{path}.kgml")
				ortho_groups_by_path[path] = 0
				for entry in tree.getroot():
					if entry.tag == 'entry':
						if entry.attrib['type'] == 'ortholog':
							ortholog = [x.replace("ko:","") for x in entry.attrib['name'].split(" ")]
							ortho_groups_by_path[path] += 1
							for graphic in entry:
								x = 0
								y = 0
								w = 0
								h = 0
								if graphic.attrib['type'] == 'rectangle':
									x = graphic.attrib['x']
									y = graphic.attrib['y']
									w = graphic.attrib['width']
									h = graphic.attrib['height']
								for ortho in ortholog:
									LOCS.write(f"{ortho}")
									LOCS.write(f"\tx:{x}")
									LOCS.write(f"\ty:{y}")
									LOCS.write(f"\tw:{w}")
									LOCS.write(f"\th:{h}")
									LOCS.write(f"\n")
			LOCS.close()


		if isfile(f"{kegg_dir}/details/ko{path}.details"):

			DET = open(f"{kegg_dir}/details/ko{path}.details",'r')

			ortho_by_path[path] = 0
			record = False
			pathway = ""
			for line in DET:
				line = line.rstrip()
				if line != "":
					if line[0:4] == "NAME":
						pathway = line[4:].strip()
						pathways[path] = pathway
					elif line[0:9] == "ORTHOLOGY":
						ortho = line[9:].strip().split(" ")[0]
						ortho_by_path[path] += 1
						if ortho not in path_by_ortho.keys(): path_by_ortho[ortho] = []
						if path not in path_by_ortho[ortho]: path_by_ortho[ortho].append(path)
						record = True
					elif line[0] == " " and record == True:
						ortho = line.strip().split(" ")[0]
						if ortho not in path_by_ortho.keys(): path_by_ortho[ortho] = []
						if path not in path_by_ortho[ortho]: path_by_ortho[ortho].append(path)
						ortho_by_path[path] += 1
					else:
						record = False

			DET.close()

	OUT = open(f"{kegg_dir}/etc/ortho_pathways.tsv",'w')
	for ortho in sorted(path_by_ortho.keys()):
		line = ";".join([f"{pathways[path]}:{ortho_by_path[path]}:{ortho_groups_by_path[path]}" for path in sorted(path_by_ortho[ortho],key = lambda x: pathways[x])])
		OUT.write(f"{ortho}\t{line}\n")
	OUT.close()
	

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()
	GetOptions.add_argument("-k","--kegg_dir",required=True)

	args = GetOptions.parse_known_args()[0]
	kegg_dir = args.kegg_dir

	kgml_parser(kegg_dir=kegg_dir)