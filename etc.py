#!/usr/bin/env	python3

def seconds_to_readable(seconds):

	m = 60
	h = 60
	d = 24

	output = ""

	days = int(seconds/(d*h*m))
	seconds = seconds%(d*h*m)

	# if days > 0:
	# 	output += f"{days}d "

	hours = int(seconds/(h*m))
	seconds = seconds%(h*m)

	# if hours > 0:
	# 	output += f"{hours}h "

	mins = int(seconds/m)
	seconds = int(seconds%m)

	# if mins > 0:
	# 	output += f"{mins}m "

	# output += f"{seconds}s"

	return f"{days:03}d {hours:02}h {mins:02}m {seconds:02}s"