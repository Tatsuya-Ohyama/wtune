#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
水分子調整プログラム
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import numpy as np
from scipy.spatial import distance as scdi
from basic_func import check_exist, check_overwrite, get_file_length
from tqdm import tqdm
import tempfile
import os

from molecule_topology import MoleculeTopology


# =============== functions =============== #
def get_shortest_distance(coord_solute, coord_solvent):
	"""
	最短距離を算出する関数
	"""
	return scdi.cdist(coord_solute, coord_solvent).min()


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
	group_strip.add_argument("-S", "--separate", dest = "separate_mode", metavar = "SEPARATE_MODE", choices = ["atom", "residue", "molecule"], default = "residue", help = "remove solvent molecule by atom, residue or molecule (Default: Residue)")
	group_strip.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly (Default: False)")

	args = parser.parse_args()

	if args.distance is None and args.number is None and args.view == False:
		sys.stderr.write("ERROR: distance or number does not specified\n")
		sys.exit(1)

	check_exist(args.input, 2)

	if args.view == True:
		# 内容表示
		with open(args.input, "r") as f_obj:
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

			print(" Input file     : %s" % args.input)
			print(" Water molecules: %5d" % w_count)
			print(" Other residues : %5d" % r_count)
			print(" All atoms      : %5d" % a_count)

	else:
		# ファイルの読み込み
		check_exist(args.input, 2)
		sys.stderr.write("Loading PDB file as topology information.")
		sys.stderr.flush()

		structure = MoleculeTopology(args.input)

		# マスクの設定
		list_solvent_index = structure.set_mask(args.mask_solvent).get_info("atom", "atom_index")
		ambermask_solute = None
		if args.mask_solute is None:
			# 溶質のマスクが指定されていない場合
			ambermask_solute = "!" + args.mask_solvent
		else:
			# 溶質のマスクが指定されている場合
			ambermask_solute = args.mask_solute

		# 溶質分子の座標取得
		list_solute_index = structure.set_mask(ambermask_solute).get_info("atom", "atom_index")
		coord_solute = np.array([[structure.molecule.atoms[atom_index].xx, structure.molecule.atoms[atom_index].xy, structure.molecule.atoms[atom_index].xz] for atom_index in list_solute_index])

		# 溶質の原子インデックスリストを作成
		list_remain_solute_idx = list(set(structure.set_mask("*").get_info("atom", "atom_index")) - set(list_solvent_index))

		sys.stderr.write("\033[2K\033[G")
		sys.stderr.flush()

		flag_first = True
		residue_info = ""
		coord_solvent = []
		list_info_distance = []	# [[residue_name, residue_number, distance, [atom_number, ...]], ...]
		for solvent_idx in tqdm(list_solvent_index, desc = "Calculate distance", ascii = True, leave = False):
			if args.separate_mode == "atom":
				# 原子レベルでのチェック
				if flag_first == False:
					# 2 回目移行は距離を計算
					list_info_distance[-1][1] = get_shortest_distance(coord_solute, np.array(coord_solvent))
				else:
					# 最初は登録するだけ
					flag_first = False

				list_info_distance.append([[], 0.0])
				coord_solvent = []

			else:
				# 残基レベルでのチェック
				res_name = structure.molecule.atoms[solvent_idx].residue.name
				res_idx = structure.molecule.atoms[solvent_idx].residue.number
				if residue_info != "{0}.{1}".format(res_name, res_idx):
					# 残基名が異なる場合
					if flag_first == False:
						# 2 回目移行は距離を計算
						list_info_distance[-1][1] = get_shortest_distance(coord_solute, np.array(coord_solvent))
					else:
						# 最初は登録するだけ
						flag_first = False

					residue_info = "{0}.{1}".format(res_name, res_idx)
					list_info_distance.append([[], 0.0])
					coord_solvent = []

			coord_solvent.append([structure.molecule.atoms[solvent_idx].xx, structure.molecule.atoms[solvent_idx].xy, structure.molecule.atoms[solvent_idx].xz])
			list_info_distance[-1][0].append(solvent_idx)

		if len(coord_solvent) != 0:
			list_info_distance[-1][1] = get_shortest_distance(coord_solute, np.array(coord_solvent))

	list_remain_solvent_idx = []
	if args.distance is not None:
		# 距離モードの場合
		for info_distance in list_info_distance:
			if info_distance[1] <= args.distance:
				list_remain_solvent_idx.extend(info_distance[0])

	elif args.number is not None:
		# 数モードの場合
		for info_distance in sorted(list_info_distance, key = lambda x : x[1]):
			list_remain_solvent_idx.extend(info_distance[0])

	if args.flag_overwrite == False:
		check_overwrite(args.output)

	atom_idx = 0
	flag_ter = False
	temp_name = ""
	max_line = get_file_length(args.input)
	flag_remain = False
	remain_mol = []
	cnt_write = 0
	with tempfile.NamedTemporaryFile(mode = "w", dir = ".", prefix = ".wtune_", delete = False) as obj_output:
		temp_name = obj_output.name
		with open(args.input, "r") as obj_input:
			for line in tqdm(obj_input, desc = "Output", ascii = True, leave = False, total = max_line):
				# PDB の読み込み
				if line.startswith("ATOM") or line.startswith("HETATM"):
					# ATOM レコードの場合
					if atom_idx in list_remain_solute_idx:
						# 溶質の書き出し
						list_remain_solute_idx.remove(atom_idx)
						obj_output.write(line)
						flag_ter = False

					elif atom_idx in list_remain_solvent_idx:
						# 残す溶媒リストの atom_id と一致する場合、書き込み
						cnt_write += 1
						list_remain_solvent_idx.remove(atom_idx)
						if args.number is not None and args.number < cnt_write:
							# 数モードの場合、上限値を越えたらリストを削除し、書き込まない
							list_remain_solvent_idx = []
							continue

						if args.separate_mode == "molecule":
							flag_remain = True
							if len(remain_mol) != 0:
								for atom_info in remain_mol:
									obj_output.write(atom_info)
								remain_mol = []
						obj_output.write(line)
						flag_ter = False

					elif args.separate_mode == "molecule":
						if flag_remain:
							obj_output.write(line)
							flag_ter = False
						else:
							remain_mol.append(line)

					atom_idx += 1

				elif line.startswith("TER"):
					# TER レコードの場合
					flag_remain = False
					remain_mol = []
					if flag_ter == False:
						obj_output.write("TER\n")
						flag_ter = True

				elif line.startswith("END") and not line.startswith("ENDMDL"):
					obj_output.write("END\n")

	os.rename(temp_name, args.output)
