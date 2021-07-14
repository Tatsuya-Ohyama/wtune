# wtune

## 概要
水分子を溶質からの特定の距離、および特定の数で切り出すプログラム


## 使用方法
```sh
$ wtune.py [-h] [-V] -i INPUT_FILE -o OUTPUT_FILE (-d DISTANCE | -n NUMBER) [-ms MASK_SOLUTE] [-mv MASK_SOLVENT] [-S SEPARATE_MODE] [-O]
```

* optional arguments:
	* `-h`, `--help`
		: ヘルプメッセージを表示して終了する。

* View mode:
	* `-V`, `--view`
		: 全原子と水分子数を表示して終了する。

* Strip mode:
	* `-i INPUT_FILE`, `--input INPUT_FILE`
		: .pdb ファイル (Input)
	* `-o OUTPUT_FILE`, `--output OUTPUT_FILE`
		: .pdb ファイル (Output)
	* `-d DISTANCE`, `--distance DISTANCE`
		: 溶質からの距離
	* `-n NUMBER`, `--number NUMBER`
		: 水分子数
	* `-ms MASK_SOLUTE`, `--mask_solute MASK_SOLUTE`
		: 溶質分子の Ambermask
	* `-mv MASK_SOLVENT`, `--mask_solvent MASK_SOLVENT`
		: 溶媒分子の Ambermask
	* `-S SEPARATE_MODE`, `--separate SEPARATE_MODE`
		: 溶媒分子を原子、残基、分子レベルで削除する (Default: Residue)
	* `-O`
		: プロンプトを出さずに上書きする (Default: False)。


## 動作要件
* Python3
	* numpy
	* scipy
	* tqdm


## License
This software is released under the MIT License, see LICENSE.


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 4.7 (2021-07-14)
* 公開した。
* PEP8 スタイルにした。
