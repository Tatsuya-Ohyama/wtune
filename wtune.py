#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
	水分子調整プログラム
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import parmed
import numpy as np
from py_module_basic import basic

# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	group_view = parser.add_argument_group("View mode")
	group_view.add_argument("-V", "--view", dest = "view", action = "store_true", help = "display the number of water molecules and all atom")

	group_strip = parser.add_argument_group("Strip mode")
	group_strip.add_argument("-i", "--input", dest = "input", required = True, help = "pdb file (input)")
	group_strip.add_argument("-o", "--output", dest = "output", required = True, help = "pdb file (output)")
	group_strip_mode = group_strip.add_mutually_exclusive_group(required = True)
	group_strip_mode.add_argument("-d", "--distance", dest = "distance", metavar = "DISTANCE", type = float, help = "specify the distance from solute")
	group_strip_mode.add_argument("-n", "--number", dest = "number", metavar = "NUMBER", type = int, help = "specify the number of water molecules")
	group_strip.add_argument("-m", "--mask", dest = "mask", default = ":SOL,WAT,HOH", help = "ambermask for solvent molecules")
	group_strip.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly (Default: False)")

	args = parser.parse_args()

	if args.distance is None and args.number is None and args.view == False:
		sys.stderr.write("ERROR: distance or number does not specified\n")
		sys.exit(1)

	basic.check_exist(args.input, 2)

	re_record_atom = re.compile(r"^((ATOM)|(HETATM))")

	if args.view == True:
		# 内容表示
		with open(args.input, "r") as f_obj:
			w_count = 0
			r_count = 0
			a_count = 0

			residue = ""
			resnum = ""
			for line in f_obj:
				if re.search(re_record_atom, line):
					a_count += 1
					if line[17:20] + line[22:26] != residue + resnum:
						residue = line[17:20]
						resnum = line[22:26]
						r_count += 1
						if line[17:20] in ["WAT", "HOH", "SOL"]:
							w_count += 1

			print(" Input file     : %s" % args.input)
			print(" Water molecules: %5d" % w_count)
			print(" Other residues : %5d" % r_count)
			print(" All atoms      : %5d" % a_count)

	else:
		# ファイルの読み込み
		basic.check_exist(args.input, 2)
		structure = parmed.load_file(args.input)
		ambermask_solvent = parmed.amber.AmberMask(structure, args.mask)
		ambermask_solute = parmed.amber.AmberMask(structure, "!" + args.mask)

		# 溶質分子の座標取得
		coord_solute = np.array([[structure.atoms[idx].xx, structure.atoms[idx].xy, structure.atoms[idx].xz] for idx in ambermask_solute.Selected()])

		# 溶質分子の残基番号取得
		remain_residue = []
		for atom_solute_idx in ambermask_solute.Selected():
			if structure.atoms[atom_solute_idx].residue.idx not in remain_residue:
				remain_residue.append(structure.atoms[atom_solute_idx].residue.idx)

		for atom_solvent_idx in ambermask_solvent.Selected():
			# 溶媒分子をチェック
			if structure.atoms[atom_solvent_idx].residue.idx not in remain_residue:
				# 溶質分子の座標をリファレンスにする
				coord_solvent = np.array([structure.atoms[atom_solvent_idx].xx, structure.atoms[atom_solvent_idx].xy, structure.atoms[atom_solvent_idx].xz])

				# 距離行列を作成
				distance_list = np.sqrt(np.sum((coord_solvent - coord_solute) ** 2, axis = 1))
				if len(np.where(distance_list <= args.distance)[0]) != 0:
					# 指定距離内に分子が存在する場合
					remain_residue.append(structure.atoms[atom_solvent_idx].residue.idx)

		strip_mask = parmed.amber.AmberMask(structure, "!:" + ",".join([str(x + 1) for x in remain_residue]))
		structure.strip(strip_mask)
		if args.flag_overwrite == False:
			basic.check_overwrite(args.output)
		structure.write_pdb(args.output, renumber = True)
