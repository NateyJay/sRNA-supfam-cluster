import sys
import os
import subprocess
import Levenshtein
import argparse
import operator
import pickle

from tqdm import tqdm

# processing inputs from command line

parser = argparse.ArgumentParser()
parser.add_argument('-file', nargs='?',required=True)
parser.add_argument('-cutoff', nargs='+', required=True)
args = parser.parse_args()

file = args.file
cutoff_list = map(int,args.cutoff)


if file[-3:] != ".fa":
	if ".fasta" not in file:
		print "FATAL: input file must be in fasta format (.fa or .fasta)."
		sys.exit()

if max(cutoff_list) > 15:
	print "FATAL: distance of 15 is the maximum accepted by the program. sRNAs more distant than even 10 are likely unrelated."
	sys.exit()

# this function defines the modified Hamming distance. It basically just tries all substrings (length = the shorter of the two strings) with the standard Hamming algorithm.
# max_cutoff is important, because larger cutoffs mean that you have to test more overhangs to find their distance. This function only finds distance up to the maximum user specified distance.

def superhamming(inpA, inpB, max_cutoff):

	if len(inpA) > len(inpB):
		A = inpB
		B = inpA
	else:
		A = inpA
		B = inpB

	len_diff = len(B) - len(A)

	best_ham = 100
	bestA = ""
	bestB = ""

	top_buffer = len_diff + max_cutoff
	bot_buffer = max_cutoff

	for i in range(top_buffer + 1):
		ii = top_buffer - i
		top_string = "." * i + A + ii * "."
		# print ""
		# print top_string

		for j in range(bot_buffer + 1):
			jj = bot_buffer - j
			bot_string = "." * j + B + jj * "."
			# print bot_string, Levenshtein.hamming(top_string, bot_string)
			ham = Levenshtein.hamming(top_string, bot_string)
			if ham < best_ham:
				best_ham = ham
				bestA = top_string
				bestB = bot_string

	return(best_ham)

clusters_by_library = {}
nodes_by_library = {}
clusters = []


# reading in fasta file

with open(file, "r") as f:
	lines = f.readlines()

for line in lines:


	line = line.strip()
	if line[0] == ">":
		header = line[1:]
	else:
		seq = line

		entry = (header, seq)

		clusters.append(entry)


# making a dictionary of relationships under the specified distance for each tested cutoff.

relationships = {}
for cutoff in cutoff_list:
	relationships[cutoff] = {}

max_cutoff = max(cutoff_list)



print '  building all by all distance matrix'
pbar = tqdm(total=len(clusters))
for i, i_entry in enumerate(clusters):
	pbar.update(1)
	for dist_cutoff in cutoff_list:
		relationships[dist_cutoff][i] = [i]


	for j, j_entry in enumerate(clusters):

		if i != j:

			i_read_name, i_seq = i_entry
			j_read_name, j_seq = j_entry

			ham = superhamming(i_seq, j_seq, max_cutoff)

			for dist_cutoff in cutoff_list:
				if ham <= dist_cutoff:
					relationships[dist_cutoff][i].append(j)

pbar.close()


# for each of those cutoffs, we find the superfamily structure

for dist_cutoff in cutoff_list:

	# first sorting the reads by the number of connections (most connected are considered first)

	print "for dist_cutoff = {}".format(dist_cutoff)
	print "   Sorting by depth..."
	keys = []
	isgrouped = {}
	for key in relationships[dist_cutoff]:
		keys.append((key,len(relationships[dist_cutoff][key])))
		isgrouped[key] = False

	keys.sort(key=operator.itemgetter(1), reverse=True)


	# Next it processes each read looking for positively grouped reads and grouping them by superfamily name. This is recursive/inclusive, meaning if a connection path can connect an sRNA to a superfamily exists, it is then grouped in that superfamily.

	print "   Building superfamily structure..."
	superfamilies = []
	for entry in keys:
		key, depth = entry

		if isgrouped[key] == False:

			if depth > 1:

				found_keys = relationships[dist_cutoff][key]

				key_buffer = list(found_keys[1:])

				while True:
					if key_buffer == []:
						break
					new_key = key_buffer.pop()

					new_found_keys = relationships[dist_cutoff][new_key]
					if len(new_found_keys) > 1:
						for y in new_found_keys[1:]:
							if y not in found_keys:
								key_buffer.append(y)
						found_keys = found_keys + new_found_keys[1:]
			else:

				found_keys = relationships[dist_cutoff][key]

			found_keys = list(set(found_keys))
			for k in found_keys:
				isgrouped[k] = True
			superfamilies.append(found_keys)


	# once the structure of the superfamilies is found, the program writes these to a file

	print "   {} <- keys".format(len(keys))
	print "   {} <- superfamilies".format(len(superfamilies))


	fam_counter = 0
	with open("01out-d{}.superfamilies.txt".format(dist_cutoff), "w") as f:
		to_print = "cluster\tsequence\tsuperfamily"

		print >> f, to_print
		
		for fam in superfamilies:
			fam_counter += 1

			fam_name = "SupFam_{}".format(fam_counter)

			entries = []
			species_contained = []
			for i in fam:
				cluster = clusters[i]

				cluster_name = cluster[0]
				cluster_seq = cluster[1]

				to_print = [cluster_name, cluster_seq, fam_name]
	
				to_print = "\t".join(to_print)

				print >> f, to_print





