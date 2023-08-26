#!/usr/bin/env	python3

def grab_KEGG_maps(outdir="KEGG_MAP_files"):

	from os.path import isdir, isfile
	from os import makedirs
	from subprocess import run
	from time import time
	from progress_bar import ProgressBar


	if not isdir(f"{outdir}/maps"):
		makedirs(f"{outdir}/maps",mode=0o755)

	if not isdir(f"{outdir}/kgmls"):
		makedirs(f"{outdir}/kgmls",mode=0o755)

	if not isdir(f"{outdir}/kgmls/ko"):
		makedirs(f"{outdir}/kgmls/ko",mode=0o755)

	if not isdir(f"{outdir}/kgmls/ec"):
		makedirs(f"{outdir}/kgmls/ec",mode=0o755)

	if not isdir(f"{outdir}/etc"):
		makedirs(f"{outdir}/etc",mode=0o755)

	if not isdir(f"{outdir}/details"):
		makedirs(f"{outdir}/details",mode=0o755)


	mapps = {}
	if isfile(f"{outdir}/etc/pathways.list"):
		PATHWAYS = open(f"{outdir}/etc/pathways.tsv",'r')
		for line in PATHWAYS:
			line = line.strip()
			mapps[line] = True
		PATHWAYS.close()
	else:
		OUT = open(f"{outdir}/etc/pathways.list",'w')
		list = run(["curl","https://rest.kegg.jp/link/module/pathway"],capture_output=True,text=True).stdout
		for line in list.split("\n"):
			if line != "":
				mapp = line.split("\t")[0].split(":")[-1]
				mapp  = mapp.replace("map","")
				if mapp not in mapps.keys():
					mapps[mapp] = True
					OUT.write(f"{mapp}\n")
		OUT.close()


	pg = ProgressBar(len(mapps.keys()))

	for mapp in mapps.keys():

		if not isfile(f"{outdir}/details/ko{mapp}.details"):
			run(["curl",f"https://rest.kegg.jp/get/ko{mapp}","-o",f"{outdir}/details/ko{mapp}.details"],capture_output=True)
		
		if not isfile(f"{outdir}/kgmls/ko/ko{mapp}.kgml"):
			run(["curl",f"https://rest.kegg.jp/get/ko{mapp}/kgml","-o",f"{outdir}/kgmls/ko/ko{mapp}.kgml"],capture_output=True)

		if not isfile(f"{outdir}/kgmls/ec/ec{mapp}.kgml"):
			run(["curl",f"https://rest.kegg.jp/get/ec{mapp}/kgml","-o",f"{outdir}/kgmls/ec/ec{mapp}.kgml"],capture_output=True)

		if not isfile(f"{outdir}/maps/ko{mapp}.png"):
			run(["wget",f"https://rest.kegg.jp/get/map{mapp}/image","-O",f"{outdir}/maps/ko{mapp}.png"],capture_output=True)

		pg.update_progress_bar()


if __name__ == "__main__":

	from argparse import ArgumentParser
	
	GetOptions = ArgumentParser()

	GetOptions.add_argument("-o","--outdir",default="KEGG_MAP_files")

	args = GetOptions.parse_known_args()[0]

	outdir = args.outdir

	grab_KEGG_maps(outdir)