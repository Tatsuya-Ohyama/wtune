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
from scipy.spatial import distance
from basic_func import check_exist, check_overwrite, summarized_range
from tqdm import tqdm


# =============== variables =============== #
re_atom = re.compile(r'^(?:(?:HETATM)|(?:ATOM))')
re_ter = re.compile(r'^TER')
re_end = re.compile(r'^END')


# =============== functions =============== #
def calc_matrix_distance(coord_solute, coord_solvent):
	sol = coord_solute.shape[0]
	solv = coord_solvent.shape[0]
	al = sol + solv

	matrix_distance = distance.pdist(np.vstack((coord_solute, coord_solvent)))

	target_list = []
	base = 0
	for i in range(sol):
		base += sol - 1 - i
		target_list.extend(list(range(base, base + solv)))
		base += solv

	return matrix_distance[np.array(target_list)].reshape(sol, solv)


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
	group_strip.add_argument("-ms", "--mask_solute", dest = "mask_solute", help = "ambermask for solvent molecules")
	group_strip.add_argument("-mv", "--mask_solvent", dest = "mask_solvent", default = ":SOL,WAT,HOH", help = "ambermask for solvent molecules")
	group_strip.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly (Default: False)")

	args = parser.parse_args()

	if args.distance is None and args.number is None and args.view == False:
		sys.stderr.write("ERROR: distance or number does not specified\n")
		sys.exit(1)

	check_exist(args.input, 2)

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
		check_exist(args.input, 2)
		structure = parmed.load_file(args.input)

		# マスクの設定
		ambermask_solvent = parmed.amber.AmberMask(structure, args.mask_solvent)
		ambermask_solute = None
		if args.mask_solute is None:
			# 溶質のマスクが指定されていない場合
			ambermask_solute = parmed.amber.AmberMask(structure, "!" + args.mask_solvent)
		else:
			# 溶質のマスクが指定されている場合
			ambermask_solute = parmed.amber.AmberMask(structure, args.mask_solute)

		# 残すマスクを作成
		remain_idx = [structure.atoms[idx].idx for idx in ambermask_solute.Selected()]

		# 溶質分子の座標取得
		coord_solute = np.array([[structure.atoms[idx].xx, structure.atoms[idx].xy, structure.atoms[idx].xz] for idx in ambermask_solute.Selected()])

		# 溶媒情報取得
		info_solvent = [[structure.atoms[idx].idx, structure.atoms[idx].residue.name, structure.atoms[idx].residue.number, structure.atoms[idx].xx, structure.atoms[idx].xy, structure.atoms[idx].xz] for idx in ambermask_solvent.Selected()]

		flag_first = True
		residue_info = ""
		atom_list = []
		coord_solvent = []
		remain_list = []	# [[residue_name, residue_number, distance, [atom_number, ...]], ...]
		for atom_record in tqdm(info_solvent, desc = "Calculate distance", ascii = True):
			if residue_info != "{0[1]}.{0[2]}".format(atom_record):
				# 残基名が異なる場合
				if flag_first == False:
					# 2 回目移行は距離を計算
					matrix_distance = calc_matrix_distance(coord_solute, np.array(coord_solvent))
					remain_list[-1][2] = matrix_distance.min()
				else:
					# 最初は登録するだけ
					flag_first = False

				residue_info = "{0[1]}.{0[2]}".format(atom_record)
				remain_list.append([atom_record[1], atom_record[2], 0.0, []])
				atom_list = []
				coord_solvent = []
			coord_solvent.append(atom_record[3:])
			remain_list[-1][3].append(atom_record[0])

		if len(coord_solvent) != 0:
			matrix_distance = calc_matrix_distance(coord_solute, np.array(coord_solvent))
			remain_list[-1][2] = matrix_distance.min()

	if args.distance is not None:
		# 距離モードの場合
		for remain_info in remain_list:
			if remain_info[2] <= args.distance:
				remain_idx.extend(remain_info[3])

	elif args.number is not None:
		# 数モードの場合
		cnt = 0
		for remain_info in sorted(remain_list, key = lambda x : x[2]):
			remain_idx.extend(remain_info[3])
			cnt += 1
			if args.number < cnt:
				break

	if args.flag_overwrite == False:
		check_overwrite(args.output)

	atom_idx = 0
	flag_ter = False
	with open(args.output, "w") as obj_output:
		with open(args.input, "r") as obj_input:
			for line in obj_input:
				# PDB の読み込み
				if re_atom.search(line):
					# ATOM レコードの場合
					if atom_idx in remain_idx:
						# 残すリストの atom_id と一致する場合、書き込み
						obj_output.write(line)
						flag_ter = False
					atom_idx += 1

				elif re_ter.search(line):
					# TER レコードの場合
					if flag_ter == False:
						obj_output.write("TER\n")
						flag_ter = True

				elif re_end.search(line):
					obj_output.write("END\n")
