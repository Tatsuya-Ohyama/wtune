#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
file function module
"""

import sys
import os
import re
import chardet



# =============== get_file_length =============== #
def get_file_length(file_path):
	"""
	Function to get number of lines in file

	Args:
		file_path (str): file path

	Returns:
		int: number of lines in file
	"""
	return sum(1 for line in open(file_path))


# =============== get_charcode =============== #
def get_charcode(input_file):
	"""
	Function to get character code of file

	Args:
		input_file (str): file path

	Returns:
		str: character code
	"""
	with open(input_file, "rb") as obj_input:
		return chardet.detect(obj_input.read())["encoding"]


# =============== get_file_list =============== #
def get_file_list(root_dir=".", max_depth=None, file_type="fd", regexp=".*"):
	"""
	function to get list of files

	Args:
		root_dir (str, optional): root directory. (Defaults: ".")
		max_depth (int, optional): depth of directory tree (Defaults: None)
		file_type (str, optional): type of files ("f": file / "d": directory / fd or df: all) (Defaults: "fd")
		regexp (str, optional): regexp for getting file. (Defaults: ".*")

	Returns:
		list: file list
	"""
	regexp = re.compile(regexp)

	if max_depth is not None:
		max_depth = int(max_depth)

	file_list = []
	for root, dirs, files in os.walk(root_dir):
		if max_depth is not None:
			tmp_path = root.replace(root_dir, "", 1)
			depth = len(tmp_path.split("/"))
			if max_depth < depth:
				continue

		if "f" in file_type:
			file_list.extend([os.path.join(root, x) for x in files if regexp.search(x)])

		if "d" in file_type:
			file_list.extend([os.path.join(root, x) for x in dirs if regexp.search(x)])

	return sorted(file_list)


# =============== normalize_filename =============== #
def convert_valid_filename(filename, replace_str="_"):
	"""
	Function to replace strings that cannot be used in file names on Windows and MacOS

	Args:
		filename (str): file name
		replace_str (str, optional): replaced char (Default: "_")

	Returns:
		str: replaced file name
	"""
	RE_PROHIBIT_CHAR = re.compile(r"[\\\/:,;\*\?\"<>\|]")
	return RE_PROHIBIT_CHAR.sub(replace_str, filename)
