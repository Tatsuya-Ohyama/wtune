#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Program to adjust water molecules
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import numpy as np
from scipy.spatial import distance as scdi
from tqdm import tqdm
import tempfile
import os

from mods.func_prompt_io import check_exist, check_overwrite
from mods.func_file import get_file_length
from mods.molecule_topology import MoleculeTopology



# =============== functions =============== #
def get_shortest_distance(coord_solute, coord_solvent):
	"""
	Function to get shortest distance

	Args:
		coord_solute (ndarray): coordinates for solute molecules
		coord_solvent (ndarray): coordinates for solvent molecules

	Returns:
		float
	"""
	return scdi.cdist(coord_solute, coord_solvent).min()



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	group_view = parser.add_argument_group("View mode")
	group_view.add_argument("-V", "--view", dest="VIEW", action="store_true", help="display the number of water molecules and all atom")

	group_strip = parser.add_argument_group("Strip mode")
	group_strip.add_argument("-i", "--input", dest="INPUT_FILE", required=True, help="pdb file (input)")
	group_strip.add_argument("-o", "--output", dest="OUTPUT_FILE", required=True, help="pdb file (output)")
	group_strip_mode = group_strip.add_mutually_exclusive_group(required=True)
	group_strip_mode.add_argument("-d", "--distance", dest="DISTANCE", metavar="DISTANCE", type=float, help="specify the distance from solute")
	group_strip_mode.add_argument("-n", "--number", dest="NUMBER", metavar="NUMBER", type=int, help="specify the number of water molecules")
	group_strip.add_argument("-ms", "--mask_solute", dest="MASK_SOLUTE", help="ambermask for solute molecules")
	group_strip.add_argument("-mv", "--mask_solvent", dest="MASK_SOLVENT", default=":SOL,WAT,HOH", help="ambermask for solvent molecules")
	group_strip.add_argument("-S", "--separate", dest="SEPARATE_MODE", metavar="SEPARATE_MODE", choices=["atom", "residue", "molecule"], default="residue", help="remove solvent molecule by atom, residue or molecule (Default: Residue)")
	group_strip.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly (Default: False)")

	args = parser.parse_args()

	if args.DISTANCE is None and args.NUMBER is None and args.VIEW == False:
		sys.stderr.write("ERROR: distance or number does not specified\n")
		sys.exit(1)

	check_exist(args.INPUT_FILE, 2)

	if args.VIEW == True:
		# display information
		with open(args.INPUT_FILE, "r") as f_obj:
			w_count = 0
			r_count = 0
			a_count = 0

			residue = ""
			resnum = ""
			for line in f_obj:
				if line.startswith("ATOM") or line.startswith("HETATM"):
					a_count += 1
					if line[17:20] + line[22:26] != residue + resnum:
						residue = line[17:20]
						resnum = line[22:26]
						r_count += 1
						if line[17:20] in ["WAT", "HOH", "SOL"]:
							w_count += 1

			print(" Input file     : %s" % args.INPUT_FILE)
			print(" Water molecules: %5d" % w_count)
			print(" Other residues : %5d" % r_count)
			print(" All atoms      : %5d" % a_count)

	else:
		# read file
		check_exist(args.INPUT_FILE, 2)
		sys.stderr.write("Loading PDB file as topology information.")
		sys.stderr.flush()

		structure = MoleculeTopology(args.INPUT_FILE)

		# set mask
		list_solvent_index = structure.set_mask(args.MASK_SOLVENT).get_info("atom", "atom_index")
		ambermask_solute = None
		if args.MASK_SOLUTE is None:
			# when Ambermask for solute is not specified
			ambermask_solute = "!" + args.MASK_SOLVENT
		else:
			# whe Ambermask for solute is specified
			ambermask_solute = args.MASK_SOLUTE

		# get coordinates for solute
		list_solute_index = structure.set_mask(ambermask_solute).get_info("atom", "atom_index")
		coord_solute = np.array([[structure.molecule.atoms[atom_index].xx, structure.molecule.atoms[atom_index].xy, structure.molecule.atoms[atom_index].xz] for atom_index in list_solute_index])

		# create index list for solute atoms
		list_remain_solute_idx = list(set(structure.set_mask("*").get_info("atom", "atom_index")) - set(list_solvent_index))

		sys.stderr.write("\033[2K\033[G")
		sys.stderr.flush()

		flag_first = True
		residue_info = ""
		coord_solvent = []
		list_info_distance = []	# [[residue_name, residue_number, distance, [atom_number, ...]], ...]
		for solvent_idx in tqdm(list_solvent_index, desc="Calculate distance", ascii=True, leave=False):
			if args.SEPARATE_MODE == "atom":
				# check atomic level
				if flag_first == False:
					# calculate distance after 2nd loop
					list_info_distance[-1][1] = get_shortest_distance(coord_solute, np.array(coord_solvent))
				else:
					# only register for 1st loop
					flag_first = False

				list_info_distance.append([[], 0.0])
				coord_solvent = []

			else:
				# check residue level
				res_name = structure.molecule.atoms[solvent_idx].residue.name
				res_idx = structure.molecule.atoms[solvent_idx].residue.number
				if residue_info != "{0}.{1}".format(res_name, res_idx):
					# when residue information differ from previous
					if flag_first == False:
						# calculate distance after 2nd loop
						list_info_distance[-1][1] = get_shortest_distance(coord_solute, np.array(coord_solvent))
					else:
						# only register for 1st loop
						flag_first = False

					residue_info = "{0}.{1}".format(res_name, res_idx)
					list_info_distance.append([[], 0.0])
					coord_solvent = []

			coord_solvent.append([structure.molecule.atoms[solvent_idx].xx, structure.molecule.atoms[solvent_idx].xy, structure.molecule.atoms[solvent_idx].xz])
			list_info_distance[-1][0].append(solvent_idx)

		if len(coord_solvent) != 0:
			list_info_distance[-1][1] = get_shortest_distance(coord_solute, np.array(coord_solvent))

	list_remain_solvent_idx = []
	if args.DISTANCE is not None:
		# distance-mode
		for info_distance in list_info_distance:
			if info_distance[1] <= args.DISTANCE:
				list_remain_solvent_idx.extend(info_distance[0])

	elif args.NUMBER is not None:
		# number-mode
		for info_distance in sorted(list_info_distance, key=lambda x : x[1]):
			list_remain_solvent_idx.extend(info_distance[0])

	if args.FLAG_OVERWRITE == False:
		check_overwrite(args.OUTPUT_FILE)

	atom_idx = 0
	flag_ter = False
	temp_name = ""
	max_line = get_file_length(args.INPUT_FILE)
	flag_remain = False
	remain_mol = []
	cnt_write = 0
	with tempfile.NamedTemporaryFile(mode="w", dir=".", prefix=".wtune_", delete=False) as obj_output:
		temp_name = obj_output.name
		with open(args.INPUT_FILE, "r") as obj_input:
			for line in tqdm(obj_input, desc="Output", ascii=True, leave=False, total=max_line):
				# read .pdb
				if line.startswith("ATOM") or line.startswith("HETATM"):
					# atom record
					if atom_idx in list_remain_solute_idx:
						# write solute information
						list_remain_solute_idx.remove(atom_idx)
						obj_output.write(line)
						flag_ter = False

					elif atom_idx in list_remain_solvent_idx:
						# write information when atom_id in remain list match atom_idx
						cnt_write += 1
						list_remain_solvent_idx.remove(atom_idx)
						if args.NUMBER is not None and args.NUMBER < cnt_write:
							# for number-mode, delete list and stop to write after over limit
							list_remain_solvent_idx = []
							continue

						if args.SEPARATE_MODE == "molecule":
							flag_remain = True
							if len(remain_mol) != 0:
								for atom_info in remain_mol:
									obj_output.write(atom_info)
								remain_mol = []
						obj_output.write(line)
						flag_ter = False

					elif args.SEPARATE_MODE == "molecule":
						if flag_remain:
							obj_output.write(line)
							flag_ter = False
						else:
							remain_mol.append(line)

					atom_idx += 1

				elif line.startswith("TER"):
					# TER record
					flag_remain = False
					remain_mol = []
					if flag_ter == False:
						obj_output.write("TER\n")
						flag_ter = True

				elif line.startswith("END") and not line.startswith("ENDMDL"):
					obj_output.write("END\n")

	os.rename(temp_name, args.OUTPUT_FILE)
