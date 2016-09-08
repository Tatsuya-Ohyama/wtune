#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
	水分子調整プログラム
"""

import argparse
import sys
import os
import re
import tqdm


# =============== functions =============== #
# ファイルの確認
def check_file(file):
	if not os.path.exists(file):
		sys.stderr.write("ERROR: No such file (%s)\n" % file)
		sys.exit(1)


# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument("-d", "--distance", help = "specify the distance from solute", metavar = "DISTANCE", type = float)
	parser.add_argument("-n", "--number", help = "specify the number of water molecules", metavar = "NUMBER", type = int)
	parser.add_argument("-t", "--thread", help = "the number of thread for calculations", metavar = "THREAD", type = int)
	parser.add_argument("-i", "--input", help = "pdb file (input)")
	parser.add_argument("-o", "--output", help = "pdb file (output)")
	parser.add_argument("-V", "--view", help = "display the number of water molecules and all atom", action = "store_true")
	parser.add_argument("-Y", "--hydrogen", help = "calculate", action = "store_true", dest = "flag_hydrogen")

	args = parser.parse_args()

	if args.distance == None and args.number == None and args.view == False:
		sys.stderr.write("ERROR: distance or number does not specified\n")
		sys.exit(1)

	check_file(args.input)

	re_record_atom = re.compile(r"^((ATOM)|(HETATM))")

	if args.view == True:
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
						if line[17:20] == "WAT" or line[17:20] == "HOH" or line[17:20] == "SOL":
							w_count += 1

			print(" %s" % args.input)
			print(" Water molecules: %5d" % w_count)
			print(" Other residues : %5d" % r_count)
			print(" All atoms      : %5d" % a_count)


	else:
		s_coords = []
		w_coords = []
		w_infos = []

		re_number = r"[\d\s]"
		with open(args.output, "w") as fobj_output:
			with open(args.input, "r") as fobj_input:
				# 溶質と水分子に分割 (溶質とその他はファイルに書き込み)

				residue = ""
				resnum = ""
				w_log = ""
				re_ter = r"^TER"
				re_wat = r"^((WAT)|(SOL)|(HOH))"
				re_end = r"^END"

				for line in fobj_input:
					# ファイル読み込み

					if re.search(re_record_atom, line):
						# 原子レコード
						if re.search(re_wat, line[17:26]):
							# 水分子の行
							if residue == line[17:26]:
								# 前と同じ残基
								w_coords[len(w_coords) - 1].append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
								w_infos[len(w_infos) - 1].append(line)

							else:
								# 前と違う残基
								w_coords.append([[float(line[30:38]), float(line[38:46]), float(line[46:54])]])
								w_infos.append([line])
								residue = line[17:26]

						else:
							# 溶質の行
							s_coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
							fobj_output.write(line)
							w_log = line

					else:
						# その他のレコード
						if not (re.search(re_ter, w_log) and re.search(re_ter, line)) and not re.search(re_end, line):
							# TER の重複回避とEND回避
							fobj_output.write(line)
							w_log = line


			# 最短距離の算出
			for i in tqdm.tqdm(range(len(w_coords))):
				# 水分子の数だけループ
				flag_break = 0

				for j in range(len(w_coords[i])):
					# 水分子構成原子の数だけループ
					saved_distance = -1

					if args.flag_hydrogen == False:
						# flag_hydrogen が1回で終了させる
						flag_break = 1

					for k in range(len(s_coords)):
						# 溶質構成原子の数だけループ
						distance = ((w_coords[i][j][0] - s_coords[k][0]) ** 2 + (w_coords[i][j][1] - s_coords[k][1]) ** 2 + (w_coords[i][j][2] - s_coords[k][2]) ** 2) ** 0.5
						if saved_distance == -1 or distance < saved_distance:
							# 初回か、これまでの距離より短い場合
							saved_distance = distance

					# 最短ルートの処理
					if 0 < saved_distance <= args.distance:
						flag_break = 1
						for k in range(len(w_infos[i])):
							fobj_output.write(w_infos[i][k])
						fobj_output.write("TER\n")

					if flag_break == 1:
						break

			fobj_output.write("END\n")



