#!/usr/bin/env python3

from time import time, sleep

class ProgressBar():

	## Number of progress ticks
	progress_ticks = 100

	## Number of divisions to split each progress bar
	divisions = 3

	## 'Progress tick' symbols
	stages = [".","\\","|","."]

	## Current step number
	step = 0
	
	## Track the average operation time
	times = []

	## How many times have been iterated through
	current_item_number = 0

	def __init__(self,items):

		## Track total time
		self.master_start = time()

		## Track total items
		self.total_items = items

		## Number of blank spaces to buffer line
		self.buffer = len(str(items))

		## Number of items between 'progress ticks'
		self.items_per_progress_tick = int(self.total_items/self.progress_ticks)

		## Number of items between 'intermediate step ticks'
		self.items_per_division = int(self.items_per_progress_tick/self.divisions)

		self.last_time = time()

		print(self.items_per_progress_tick,self.items_per_division)


		print("\n\n\n")

	def update_progress_bar(self):

		self.current_item_number += 1

		current_progress_tick = int(self.current_item_number/self.items_per_progress_tick)

		current_division = int((self.current_item_number%self.items_per_progress_tick)/self.items_per_division)

		complete = "/"*current_progress_tick

		ticks_remaining = self.progress_ticks - current_progress_tick

		incomplete = "."*ticks_remaining

		bar = complete+ self.stages[current_division] + incomplete

		progress = self.current_item_number/self.total_items*100

		if len(self.times) == 100:
			self.times.pop(0)
		

		elapsed = time() - self.last_time
		self.times.append(elapsed)

		items_remaining = self.total_items-self.current_item_number
		time_remaining = items_remaining*((sum(self.times))/len(self.times))

		printed_time = self.seconds_to_readable(time_remaining)

		up = "\x1B[3A"

		total_elapsed = self.seconds_to_readable(time()-self.master_start)

		print(f"{up}\t{self.current_item_number: {self.buffer}}/{self.total_items}  [{bar}]  {progress: 3.2f}%")
		print(f"\t\t\t  Time Elapsed: {total_elapsed}")
		print(f"\t\t\tTime Remaining: {printed_time}")

		self.last_time = time()

	def seconds_to_readable(self,seconds):

		m = 60
		h = 60
		d = 24

		output = ""

		days = int(seconds/(d*h*m))
		seconds = seconds%(d*h*m)

		hours = int(seconds/(h*m))
		seconds = seconds%(h*m)

		mins = int(seconds/m)
		seconds = int(seconds%m)

		return f"{days:03}d {hours:02}h {mins:02}m {seconds:02}s"


if __name__ == "__main__":

	items = 10000

	progress = ProgressBar(items)

	for i in range(items):
		
		progress.update_progress_bar()

		sleep(0.1)