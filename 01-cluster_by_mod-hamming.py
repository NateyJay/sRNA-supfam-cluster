import sys
import os
import subprocess
import Levenshtein
import argparse
# import networkx as nx
# import matplotlib.pyplot as plt
# import itertools
import operator
import pickle

from tqdm import tqdm


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
# num_files = len(cluster_files)
# num_files = 7

# file = '/Volumes/Keep/+Cuscuta_evolution/NextSeq_Sep2018_sRNAseq//clustering_trimmed/04out-all.clusters.HI.txt'
with open(file, "r") as f:
	lines = f.readlines()

for line in lines:
# for line in lines[1:500]:

	line = line.strip()
	if line[0] == ">":
		header = line[1:]
	else:
		seq = line

		entry = (header, seq)

		clusters.append(entry)

# for i, file in enumerate(cluster_files):
# 	print i, file

# 	clusters_by_library[species[i]] = []
# 	nodes_by_library[i] = []

# 	with open(file, "r") as f:
# 		for line in f:
# 			line = line.strip()
# 			if line[0] == ">":
# 				read_name = line[1:]
# 			else:
# 				seq = line

# 				entry = (read_name, seq)
# 				clusters.append(entry)
# 				clusters_by_library[species[i]].append(entry)

# print len(clusters)
# sys.exit()

relationships = {}
for cutoff in cutoff_list:
	relationships[cutoff] = {}

max_cutoff = max(cutoff_list)

# relationships = {0:{}, 2:{}, 5:{}}
# for dist_cutoff in [0,2,5]:
# 	relationships[dist_cutoff] = {}
# clusters = clusters[:200]

# all_file = "00sub-allbyall.dict"

# if all_ is True:
# 	print 'all by all dict found!'
# 	with open(all_file, 'r') as d:
# 		relationships = pickle.load(d)

# else:
# 	print 'all by all dict not found...'
print '  building all by all distance matrix'
pbar = tqdm(total=len(clusters))
for i, i_entry in enumerate(clusters):
	pbar.update(1)
	for dist_cutoff in cutoff_list:
		relationships[dist_cutoff][i] = [i]
	# relationships[i] = [i]
	for j, j_entry in enumerate(clusters):

		if i != j:

			i_read_name, i_seq = i_entry
			j_read_name, j_seq = j_entry

			ham = superhamming(i_seq, j_seq, max_cutoff)

			for dist_cutoff in cutoff_list:
				if ham <= dist_cutoff:
					relationships[dist_cutoff][i].append(j)

pbar.close()
# print '  writing to file "00sub-allbyall.dict"'
# with open(all_file, 'w') as d:
# 	pickle.dump(relationships, d)

for dist_cutoff in cutoff_list:
	print "for dist_cutoff = {}".format(dist_cutoff)
	print "   Sorting by depth..."
	keys = []
	isgrouped = {}
	for key in relationships[dist_cutoff]:
		keys.append((key,len(relationships[dist_cutoff][key])))
		isgrouped[key] = False

	keys.sort(key=operator.itemgetter(1), reverse=True)
	# print


	print "   Building superfamily structure..."
	superfamilies = []
	for entry in keys:
		key, depth = entry

		if isgrouped[key] == False:

			if depth > 1:

				found_keys = relationships[dist_cutoff][key]


				key_buffer = list(found_keys[1:])
				# print "KB:",key_buffer

				while True:
					if key_buffer == []:
						break
					new_key = key_buffer.pop()
					# print "KB:", key_buffer
					# print new_key


					new_found_keys = relationships[dist_cutoff][new_key]
					if len(new_found_keys) > 1:
						for y in new_found_keys[1:]:
							if y not in found_keys:
								key_buffer.append(y)
						found_keys = found_keys + new_found_keys[1:]
			else:

				found_keys = relationships[dist_cutoff][key]

				# print "KB:",key_buffer



			found_keys = list(set(found_keys))
			for k in found_keys:
				isgrouped[k] = True
			superfamilies.append(found_keys)
		# for i in found_keys:
		# 	print i
		# 	del relationships[i]

	print "   {} <- keys".format(len(keys))
	print "   {} <- superfamilies".format(len(superfamilies))

	species = ["cca","cgr-dp","cgr-pm","cgr-mass","cpe-2015","cpe-2017","cind"]
	fam_counter = 0
	with open("01out-d{}.superfamilies.txt".format(dist_cutoff), "w") as f:
		to_print = "cluster\tsequence\tspecies\tsuperfamily\teulercode"
		for spec in species:
			to_print = "{}\t{}".format(to_print, spec)
		print >> f, to_print
		for fam in superfamilies:
			fam_counter += 1

			fam_name = "SupFam_{}".format(fam_counter)

			entries = []
			species_contained = []
			for i in fam:
				cluster = clusters[i]
				# print cluster


				# print cluster
				cluster_name = cluster[0]
				cluster_species = cluster[0].split("_")[1]


				species_contained.append(cluster_species)
				entry = (cluster_name, cluster[1], cluster_species)
				entries.append(entry)


			euler_key = []

			for spec in species:
				if spec in species_contained:
					euler_key.append("1")
				else:
					euler_key.append("0")
			species_contained = list(set(species_contained))

			species_contained.sort()

			for entry in entries:

				spec_id = "&".join(species_contained)

				# euler_key = "\t".join(map(str, euler_key))

				to_print = list(entry) + [fam_name, spec_id] + euler_key
				to_print = "\t".join(to_print).replace('cca','ccm').replace('cind','cin')

				# print to_print
				print >> f, to_print

















