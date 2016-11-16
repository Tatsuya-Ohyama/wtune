#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
	水分子調整プログラム
"""

import sys, os, re, signal
signal.signal(signal.SIGINT, signal.SIG_DFL)

import argparse
import joblib
import tempfile, shutil

# =============== functions =============== #
# ファイルの確認
def check_file(file):
	if not os.path.exists(file):
		sys.stderr.write("ERROR: No such file (%s)\n" % file)
		sys.exit(1)


# check_overwrite (上書き確認)
def check_overwrite(file):
	if os.path.exists(file):
		sys.stderr.write("WARN: %s exists. Overwrite it? (y/N): " % file)
		sys.stderr.flush()
		user = sys.stdin.readline().replace("\n", "")
		if not user == "y": #or not user == "Y":
			sys.exit(0)


# 溶質からの最短距離算出
def calc_sdistance(index):
	results = []
	flag_break = 0
	saved_distance = -1
	for j in range(len(w_coords[index])):
		# 水分子構成原子の数だけループ
		if args.flag_hydrogen == False:
			# flag_hydrogen が OFF の場合、1回で終了させる
			flag_break = 1

		for k in range(len(s_coords)):
			# 溶質構成原子の数だけループ
			distance = ((w_coords[index][j][0] - s_coords[k][0]) ** 2 + (w_coords[index][j][1] - s_coords[k][1]) ** 2 + (w_coords[index][j][2] - s_coords[k][2]) ** 2) ** 0.5
			if saved_distance == -1 or distance < saved_distance:
				# 初回か、これまでの距離より短い場合
				saved_distance = distance

		if flag_break == 1:
			break

	# 最短ルートの処理
	if 0 < saved_distance:
		results.append(saved_distance)
		results.append(index)
		for k in range(len(w_infos[index])):
			results.append(w_infos[index][k])
		results.append("TER\n")
	else:
		sys.stderr.write("ERROR: function (calc_sdistance) error\n")
		sys.exit(1)

	return results



# =============== main =============== #
if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	group = parser.add_mutually_exclusive_group(required = True)
	group.add_argument("-d", "--distance", dest = "distance", metavar = "DISTANCE", type = float, help = "specify the distance from solute")
	group.add_argument("-n", "--number", dest = "number", metavar = "NUMBER", type = int, help = "specify the number of water molecules")
	parser.add_argument("-t", "--thread", dest = "thread", metavar = "THREAD", type = int, default = 1, help = "the number of thread for calculations (Default: 1)")
	parser.add_argument("-i", "--input", dest = "input", required = True, help = "pdb file (input)")
	parser.add_argument("-o", "--output", dest = "output", required = True, help = "pdb file (output)")
	parser.add_argument("-V", "--view", dest = "view", action = "store_true", help = "display the number of water molecules and all atom")
	parser.add_argument("-O", dest = "flag_overwrite", action = "store_true", default = False, help = "overwrite forcibly (Default: False)")
	parser.add_argument("-Y", "--hydrogen", dest = "flag_hydrogen", action = "store_true", help = "calculate distance between solute and water molecules with hydrogen atoms (Default: False)")

	args = parser.parse_args()

	if args.distance == None and args.number == None and args.view == False:
		sys.stderr.write("ERROR: distance or number does not specified\n")
		sys.exit(1)

	check_file(args.input)

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
						if line[17:20] == "WAT" or line[17:20] == "HOH" or line[17:20] == "SOL":
							w_count += 1

			print(" Input file     : %s" % args.input)
			print(" Water molecules: %5d" % w_count)
			print(" Other residues : %5d" % r_count)
			print(" All atoms      : %5d" % a_count)

	else:
		# 水分子切り出し
		s_coords = []
		w_coords = []
		w_infos = []

		check_overwrite(args.output)

		re_number = r"[\d\s]"
		tempfile_name = ""
		with tempfile.NamedTemporaryFile(mode = "w", prefix = ".wtune_", delete = False) as obj_output:
			tempfile_name = obj_output.name
			with open(args.input, "r") as obj_input:
				# 溶質と水分子に分割 (溶質とその他はファイルに書き込み)

				residue = ""
				resnum = ""
				w_log = ""
				re_ter = r"^TER"
				re_wat = r"^((WAT)|(SOL)|(HOH))"
				re_end = r"^END"

				for line in obj_input:
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
							obj_output.write(line)
							w_log = line

					else:
						# その他のレコード
						if not (re.search(re_ter, w_log) and re.search(re_ter, line)) and not re.search(re_end, line):
							# TER の重複回避とEND回避
							obj_output.write(line)
							w_log = line


			# 最短距離の算出
			datas = joblib.Parallel(n_jobs = args.thread, verbose = 10)([joblib.delayed(calc_sdistance)(i) for i in range(len(w_coords))])

			if args.distance != None:
				for data in datas:
					if data[0] <= args.distance:
						for i in range(2, len(data)):
							obj_output.write(data[i])

			else:
				datas = sorted(datas, key = lambda x:x[0])
				datas = datas[0:args.number]
				datas = sorted(datas, key = lambda x:x[1])
				for data in datas:
					for i in range(2, len(data)):
						obj_output.write(data[i])

			obj_output.write("END\n")

		shutil.move(tempfile_name, args.output)
