#!/usr/bin/env python3

name = 'kgml_parser.py'
updated = '2023-08-02'
version = '0.2.0'

def kgml_parser(kgml_dir=None,path_list=None,outdir="KEGG_pathway_loc_files"):

	from os.path import isdir, isfile, basename
	from os import makedirs, listdir, stat
	import xml.etree.ElementTree as ET


	if not isdir(outdir):
		makedirs(outdir,mode=0o755)

	if not isdir(f"{outdir}/ec"):
		makedirs(f"{outdir}/ec",mode=0o755)

	if not isdir(f"{outdir}/ko"):
		makedirs(f"{outdir}/ko",mode=0o755)


	paths = {}
	PL = open(path_list,'r')
	for line in PL:
		paths[line.strip()] = True
	PL.close()


	locs = {}

	for path in paths.keys():

		if isfile(f"{kgml_dir}/ec/ec{path}.kgml"):

			LOCS = open(f"{outdir}/ec/ec{path}.locs",'w')
			if stat(f"{kgml_dir}/ec/ec{path}.kgml").st_size != 0:
			
				tree = ET.parse(f"{kgml_dir}/ec/ec{path}.kgml")

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

		if isfile(f"{kgml_dir}/ko/ko{path}.kgml"):

			LOCS = open(f"{outdir}/ko/ko{path}.locs",'w')
			if stat(f"{kgml_dir}/ko/ko{path}.kgml").st_size != 0:
				tree = ET.parse(f"{kgml_dir}/ko/ko{path}.kgml")
				for entry in tree.getroot():
					if entry.tag == 'entry':
						if entry.attrib['type'] == 'ortholog':
							ortholog = [x.replace("ko:","") for x in entry.attrib['name'].split(" ")]
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

if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()
	GetOptions.add_argument("-k","--kgml_dir",required=True)
	GetOptions.add_argument("-p","--path_list",required=True)
	GetOptions.add_argument("-o","--outdir",default="KEGG_pathway_loc_files")

	args = GetOptions.parse_known_args()[0]
	kgml_dir = args.kgml_dir
	path_list = args.path_list
	outdir = args.outdir

	kgml_parser(kgml_dir=kgml_dir,path_list=path_list,outdir=outdir)