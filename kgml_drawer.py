#!/usr/bin/env python3

name = 'kgml_drawer.py'
updated = '2023-08-12'
version = '0.5.1'

def kgml_drawer(custom_accessions=None,utk_dir=None,kegg_files_dir=None,kegg_locs_dir=None,map_aio=False,outdir="Filled_KEGG_Maps"):


	from PIL import Image, ImageDraw, ImageFont
	from os.path import isfile, isdir
	from os import makedirs


	if not isdir(outdir):
		makedirs(outdir,mode=0o755)

	if map_aio and not isdir(f"{outdir}/AIO"):
		makedirs(f"{outdir}/AIO",mode=0o755)

	ecs = {}
	kos = {}

	CONV = open(f"{utk_dir}/conversions.tsv",'r')
	for line in CONV:
		if line != "":

			line = line.strip()
			ncbi,upas,pae,ko,ec = line.split("\t")

			ecs[ncbi] = ec
			kos[ncbi] = ko

	CONV.close()


	by_sample = {}

	if not custom_accessions:
		
		ABS = open(f"{utk_dir}/accession_by_sample.tsv",'r')
		for line in ABS:
			
			line = line.strip()
			data = line.split("\t")
			ncbi = data[0]

			for sample in data[1:]:

				sample,count = sample.split(":")

				if sample not in by_sample.keys():
					by_sample[sample] = {}

				count = int(count)

				if count == 1:

					ko = kos[ncbi]
					ec = ecs[ncbi]

					by_sample[sample][ko] = True
					by_sample[sample][ec] = True
		ABS.close()

	else:

		CUS = open(f"{custom_accessions}",'r')
		samples = None
		for line in CUS:

			line = line.strip()

			if not samples:
				
				samples = line.split("\t")

			else:

				data = line.split("\t")

				for index,ncbi in enumerate(data):

					sample = samples[index]

					if ncbi != "":

						if sample not in by_sample.keys():
							by_sample[sample] = {}

						ko = kos[ncbi]
						ec = ecs[ncbi]

						by_sample[sample][ko] = True
						by_sample[sample][ec] = True

		CUS.close()


	pathways = {}
	PATHS = open(f"{kegg_files_dir}/etc/pathways.list",'r')
	for line in PATHS:
		pathways[line.strip()] = True
	PATHS.close()


	for path in pathways.keys():

		## Track all protein coordinates, and which KO/EC are assigned to them
		coords = {}
		
		## Keep track of the coordinates for each KO/EC
		rectangles = {}

		## Loading in Enzyme Commision locations
		if isfile(f"{locs}/ec/ec{path}.locs"):
			ECL = open(f"{locs}/ec/ec{path}.locs",'r')
			for line in ECL:
				line = line.strip()
				ec,x,y,w,h = line.split("\t")
				x = int(x.split(":")[-1])
				y = int(y.split(":")[-1])
				h = int(h.split(":")[-1])
				w = int(w.split(":")[-1])
				x1 = x - 0.5*w
				x2 = x + 0.5*w
				y1 = y - 0.5*h
				y2 = y + 0.5*h
				coor = f"{x1}:{y1},{x2}:{y2}"

				if ec not in rectangles.keys():
					rectangles[ec] = []
				rectangles[ec].append(coor)

				if coor not in coords.keys():
					coords[coor] = {}

				coords[coor][ec] = True

			ECL.close()

		## Loading in KEGG Ortholog locations
		if isfile(f"{locs}/ko/ko{path}.locs"):
			KOL = open(f"{locs}/ko/ko{path}.locs",'r')
			for line in KOL:
				line = line.strip()
				ko,x,y,w,h = line.split("\t")
				x = int(x.split(":")[-1])
				y = int(y.split(":")[-1])
				h = int(h.split(":")[-1])
				w = int(w.split(":")[-1])
				x1 = x - 0.5*w
				x2 = x + 0.5*w
				y1 = y - 0.5*h
				y2 = y + 0.5*h
				coor = f"{x1}:{y1},{x2}:{y2}"

				if ko not in rectangles.keys():
					rectangles[ko] = []
				rectangles[ko].append(coor)

				if coor not in coords.keys():
					coords[coor] = {}

				coords[coor][ko] = True
				
			KOL.close()

		if map_aio:

			colors = [
					(255,000,000,100), ## RED
					(000,255,000,100), ## GREEN
					(000,000,255,100), ## BLUE
					(255,255,000,100), ## CYAN
				]

			## Create a legend
			labels = [x for x in sorted(by_sample.keys())]
			sizes = [len(x) for x in labels]
			longest_label = labels[sizes.index(max(sizes))]
			font = ImageFont.truetype("LiberationMono-Regular.ttf", size=1600)
			text_width,text_height = font.getsize(longest_label)

			box_height = 17*100
			box_width = 46*100

			height = int(box_height*len(by_sample.keys())*1.5)
			width = int(box_width*1.5) + int(text_width)
			legend = Image.new('RGBA',(width,height),(255,255,255,0))
			draw = ImageDraw.Draw(legend)

			x1 = int(box_height*0.125)
			y1 = int(box_width*0.125)
			for index,sample in enumerate(by_sample.keys()):
				draw.rectangle([x1,y1,x1+box_width,y1+box_height],colors[index])
				draw.text((x1+int(box_width*1.25),y1),labels[index],(0,0,0),font=font)
				y1 += int(box_height*1.5)


			legend.save(f"{outdir}/AIO/legend.png")


			## Open the KEGG map, create a blank canvas to draw over
			image = Image.open(f"{kegg_files_dir}/maps/ko{path}.png")
			new = Image.new('RGBA',image.size,(255,255,255,0))
			draw = ImageDraw.Draw(new)

			## Iterate over each coordinate
			for coor in coords.keys():
				
				c1,c2 = coor.split(",")
				x1,y1 = c1.split(":")
				x2,y2 = c2.split(":")

				x1 = float(x1)
				y1 = float(y1)
				x2 = float(x2)
				y2 = float(y2)

				width = x2 - x1

				color = (0,0,0,100)

				total = len(by_sample.keys())

				## Check which samples has the protein present at the coord
				present = {x:True for x in by_sample.keys() for y in coords[coor].keys() if y in by_sample[x].keys()}


				## Divide the rectangle into pieces!
				for index,sam in enumerate(sorted(by_sample.keys())):

					temp_color = color

					x2 = x1 + (width/total)

					true = False
					if sam in present.keys():
						temp_color = colors[index]
						true = True
					
					draw.rectangle([x1,y1,x2,y2],fill=temp_color)

					x1 += (width/total)


			out = Image.alpha_composite(image,new)

			out.save(f"{outdir}/AIO/ko{path}.png")

		
		## Sample by Sample
		for sam in by_sample.keys():

			temp_dir = f"{outdir}/{sam}"

			if not isdir(temp_dir):
				makedirs(temp_dir,mode=0o755)

			image = Image.open(f"{kegg_files_dir}/maps/ko{path}.png")
			new = Image.new('RGBA',image.size,(255,255,255,0))
			draw = ImageDraw.Draw(new)

			for coor in coords.keys():
				
				c1,c2 = coor.split(",")
				x1,y1 = c1.split(":")
				x2,y2 = c2.split(":")

				x1 = float(x1)
				y1 = float(y1)
				x2 = float(x2)
				y2 = float(y2)

				color = (0,0,0,100)

				total = len(by_sample.keys())
		
				for coord_locator in coords[coor]:

					if coord_locator in by_sample[sam].keys():
						
						color = (255,0,255,100)

				draw.rectangle([x1,y1,x2,y2],fill=color)


			out = Image.alpha_composite(image,new)

			out.save(f"{outdir}/{sam}/ko{path}.png")


		## Complex coloring
		# for sam in by_sample.keys():

		# 	temp_dir = f"{outdir}/{sam}"

		# 	if not isdir(temp_dir):
		# 		makedirs(temp_dir,mode=0o755)

		# 	image = Image.open(f"{kegg_files_dir}/maps/ko{path}.png")
		# 	new = Image.new('RGBA',image.size,(255,255,255,0))
		# 	draw = ImageDraw.Draw(new)

		# 	for coor in coords_colors.keys():
				
		# 		c1,c2 = coor.split(",")
		# 		x1,y1 = c1.split(":")
		# 		x2,y2 = c2.split(":")

		# 		x1 = float(x1)
		# 		y1 = float(y1)
		# 		x2 = float(x2)
		# 		y2 = float(y2)

		# 		color = coords_colors[coor]

		# 		total = len(by_sample.keys())
		# 		present = {x:0 for x in by_sample.keys()}
		
		# 		for coord_locator in coords[coor]:

		# 			for sam_2 in by_sample.keys():

		# 				if coord_locator in by_sample[sam_2].keys():

		# 					present[sam_2] = 1


		# 		count = sum(present.values())

		# 		## Not present in any sample, black
		# 		if count == 0:
		# 			color = (0,0,0,100)

		# 		## Present in every sample, cyan
		# 		elif total == count:
		# 			color = (0,255,255,100)

		# 		## Present in every sample but this one, red
		# 		elif present[sam] == 0 and count == total - 1:
		# 			color = (255,0,0,100)
				
		# 		## Present in only this one, magenta
		# 		elif present[sam] == 1 and count == 1:
		# 			color = (255,255,0,100)
				
		# 		## Present in this sample, green
		# 		elif present[sam] == 1:
		# 			color = (0,255,0,100)

		# 		## Missing in this sample, orange
		# 		elif present[sam] == 0:
		# 			color = (255,102,0,100)


		# 		draw.rectangle([x1,y1,x2,y2],fill=color)


		# 	out = Image.alpha_composite(image,new)

		# 	out.save(f"{outdir}/{sam}/ko{path}.png")

		
if __name__ == "__main__":

	from argparse import ArgumentParser

	GetOptions = ArgumentParser()
	GetOptions.add_argument("-c","--custom_accessions_list",default=False)
	GetOptions.add_argument("-k","--kegg_files_dir",required=True)
	GetOptions.add_argument("-u","--utk_dir",required=True)
	GetOptions.add_argument("-l","--locs_dir",required=True)
	GetOptions.add_argument("-a","--all_in_one",default=False,action='store_true')
	GetOptions.add_argument("-o","--outdir",default="FILLED_MAPS")

	args = GetOptions.parse_known_args()[0]

	custom_accessions = args.custom_accessions_list
	kegg_dir = args.kegg_files_dir
	utk_dir = args.utk_dir
	locs = args.locs_dir
	aio = args.all_in_one
	outdir = args.outdir

	kgml_drawer(
		custom_accessions=custom_accessions,
		kegg_files_dir=kegg_dir,
		utk_dir=utk_dir,
		kegg_locs_dir=locs,
		map_aio=aio,outdir=outdir
	)