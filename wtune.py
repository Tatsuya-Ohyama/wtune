#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Program to adjust water molecules
"""

import sys, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import numpy as np
from scipy.spatial import distance as scdi
import parmed

from mods.func_prompt_io import check_exist, check_overwrite



# =============== constant =============== #
WATER_RESIDUES = ["SOL", "WAT", "HOH"]


# =============== functions =============== #
def view_structure(input_file):
	"""
	function to display structure summary

	Args:
		obj_mol (obj_molecule): parmed molecule object
	"""
	obj_mol = parmed.load_file(input_file)
	list_residue = [obj_residue.name for obj_residue in obj_mol.residues]
	print("{0}]".format(input_file))
	print(" * {0:<18}: {1}".format("Number of atoms", len(obj_mol.atoms)))
	print(" * {0:<18}: {1}".format("Number of residues", len(obj_mol.residues)))
	print(" * {0:<18}: {1}".format("Water molecules", len([res_name for res_name in list_residue if res_name in WATER_RESIDUES])))
	print(" * {0:<18}: {1}".format("Other residues", len([res_name for res_name in list_residue if res_name not in WATER_RESIDUES])))


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
	subparser = parser.add_subparsers()

	common_argument_group = argparse.ArgumentParser(add_help=False)
	common_argument_group.add_argument("-i", dest="INPUT_FILE", metavar="INPUT.pdb", required=True, help="pdb file (input)")

	subparser_view = subparser.add_parser("view", description="display structure summary", parents=[common_argument_group])
	subparser_view.set_defaults(handler="view")

	subparser_extract = subparser.add_parser("extract", description="extract solvent molecules", parents=[common_argument_group])
	subparser_extract.set_defaults(handler="extract")
	subparser_extract.add_argument("-o", "--output", dest="OUTPUT_FILE", required=True, help="pdb file (output)")
	subparser_extract_mode = subparser_extract.add_mutually_exclusive_group(required=True)
	subparser_extract_mode.add_argument("-d", "--distance", dest="DISTANCE", metavar="DISTANCE", type=float, help="distance from solute")
	subparser_extract_mode.add_argument("-n", "--number", dest="NUMBER", metavar="NUMBER", type=int, help="the number of solvent molecules")
	subparser_extract.add_argument("-ms", "--mask_solute", dest="MASK_SOLUTE", help="Ambermask for solute molecules")
	subparser_extract.add_argument("-mv", "--mask_solvent", dest="MASK_SOLVENT", default=":"+",".join(WATER_RESIDUES), help="Ambermask for solvent molecules (Default: `:{0}`)".format(",".join(WATER_RESIDUES)))
	subparser_extract.add_argument("-S", "--separate", dest="SEPARATE_MODE", metavar="SEPARATE_MODE", choices=["atom", "residue"], default="residue", help="remove solvent molecule by `atom` or `residue` (Default: residue)")
	subparser_extract.add_argument("-O", dest="FLAG_OVERWRITE", action="store_true", default=False, help="overwrite forcibly (Default: False)")

	args = parser.parse_args()

	check_exist(args.INPUT_FILE, 2)

	if args.handler == "view":
		# view mode
		view_structure(args.INPUT_FILE)

	elif args.handler == "extract":
		# extract mode
		if args.DISTANCE is None and args.NUMBER is None:
			sys.stderr.write("ERROR: distance or number does not specified\n")
			sys.exit(1)

		# load molecule
		obj_mol = parmed.load_file(args.INPUT_FILE)

		# set mask
		list_atom_number_all = [obj_atom.number for obj_atom in obj_mol.atoms]
		obj_ambermask_solvent = parmed.amber.AmberMask(obj_mol, args.MASK_SOLVENT)
		list_atom_idx_solvent = [obj_atom for obj_atom in obj_ambermask_solvent.Selected()]
		ambermask_solute = None
		if args.MASK_SOLUTE is None:
			# when Ambermask for solute is not specified
			ambermask_solute = "!({0})".format(args.MASK_SOLVENT)
		else:
			# whe Ambermask for solute is specified
			ambermask_solute = args.MASK_SOLUTE
		obj_ambermask_solute = parmed.amber.AmberMask(obj_mol, ambermask_solute)
		list_atom_idx_solute = list(obj_ambermask_solute.Selected())

		# get coordinates
		coord_solute = np.array([[obj_mol.atoms[atom_index].xx, obj_mol.atoms[atom_index].xy, obj_mol.atoms[atom_index].xz] for atom_index in list_atom_idx_solute])
		coord_solvent = np.array([[obj_mol.atoms[atom_index].xx, obj_mol.atoms[atom_index].xy, obj_mol.atoms[atom_index].xz] for atom_index in list_atom_idx_solvent])

		# calculate distance
		matrix_distance = scdi.cdist(coord_solute, coord_solvent)	# row: solute / column: solvent

		list_atom_idx_extract = []
		if args.DISTANCE is not None:
			# distance mode
			column_idx_within = np.where(matrix_distance <= args.DISTANCE)[1]
			for column_idx in column_idx_within:
				atom_idx = list_atom_idx_solvent[column_idx]
				# atom_name = obj_mol.atoms[atom_idx].name
				# atom_number = obj_mol.atoms[atom_idx].number
				list_atom_idx_extract.append(atom_idx)

			if args.SEPARATE_MODE.lower() == "residue":
				# residue-unit
				set_atom_idx_extract = set(list_atom_idx_extract)
				atom_idx_s = 0
				for obj_residue in obj_mol.residues:
					n_atom = len(obj_residue.atoms)
					set_atom_idxs = set(range(atom_idx_s, atom_idx_s + n_atom))
					cross_atom = set_atom_idxs & set_atom_idx_extract
					if len(cross_atom) != 0:
						set_atom_idx_extract |= set_atom_idxs
					atom_idx_s += n_atom
				list_atom_idx_extract = list(sorted(set_atom_idx_extract))

		elif args.NUMBER is not None:
			# number mode

			# minimum distance
			dist_min = np.sort(matrix_distance, axis=0)[0]
			columns_idxs = dist_min.argsort().tolist()

			if args.SEPARATE_MODE.lower() == "atom":
				# atom-unit
				list_atom_idx_extract = [list_atom_idx_solvent[column_idx] for column_idx in columns_idxs[:args.NUMBER]]

			else:
				# residue-unit
				n_residue = 0
				set_atom_idx_extract = set([])
				for column_idx in columns_idxs:
					atom_idx = list_atom_idx_solvent[column_idx]
					list_atom_number_extract = [obj_atom.number for obj_atom in obj_mol.atoms[atom_idx].residue.atoms]
					set_atom_idxs = set([list_atom_number_all.index(atom_number) for atom_number in list_atom_number_extract])
					if len(set_atom_idx_extract | set_atom_idxs) > len(set_atom_idx_extract):
						set_atom_idx_extract |= set_atom_idxs
						n_residue += 1
						if n_residue >= args.NUMBER:
							break

				list_atom_idx_extract = list(sorted(set_atom_idx_extract))

		remain_atom_idx = list_atom_idx_solute + list_atom_idx_extract

		# delete atoms
		list_mask = [1 for v in obj_mol.atoms]
		obj_mol.strip([0 if i in remain_atom_idx else 1 for i in range(len(list_mask))])

		if args.FLAG_OVERWRITE == False:
			check_overwrite(args.OUTPUT_FILE)

		obj_mol.write_pdb(args.OUTPUT_FILE, increase_tercount=True)
