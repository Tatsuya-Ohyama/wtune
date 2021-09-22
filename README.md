# wtune

## 概要
溶媒分子を溶質からの特定の距離、および特定の数で切り出すプログラム


## 使用方法
### View mode
全原子と水分子数を表示して終了する。

```sh
$ wtune.py view [-h] -i INPUT.pdb
```

* `-h`, `--help`
	: ヘルプメッセージを表示して終了する。
* `-i INPUT.pdb`
	: .pdb ファイル (Input)

### Extract mode
溶媒分子を距離あるいは数で切り出す。

```sh
$ wtune.py extract [-h] -i INPUT.pdb -o OUTPUT_FILE (-d DISTANCE | -n NUMBER) [-ms MASK_SOLUTE] [-mv MASK_SOLVENT] [-S SEPARATE_MODE] [-O]
```

* `-h`, `--help`
	: ヘルプメッセージを表示して終了する。
* `-i INPUT.pdb`
	: .pdb ファイル (Input)
* `-o OUTPUT_FILE`, `--output OUTPUT_FILE`
	: .pdb ファイル (Output)
* `-d DISTANCE`, `--distance DISTANCE`
	: 切り出す溶媒分子を溶質からの距離で指定
* `-n NUMBER`, `--number NUMBER`
	: 切り出す溶媒分子を溶質から近い順に分子数で指定
* `-ms MASK_SOLUTE`, `--mask_solute MASK_SOLUTE`
	: 溶質分子の Ambermask
* `-mv MASK_SOLVENT`, `--mask_solvent MASK_SOLVENT`
	: 溶媒分子の Ambermask (Default: `:SOL,WAT,HOH`)
* `-S SEPARATE_MODE`, `--separate SEPARATE_MODE`
	: 切り出す対象モード (`atom` or `residue`) (Default: `residue`)
* `-O`
	: プロンプトを出さずに上書き


## 動作要件
* Python3
	* numpy
	* scipy
	* parmed


## License
This software is released under the MIT License, see LICENSE.


## Authors
* Tatsuya Ohyama


## ChangeLog
### Ver. 5.1 (2021-09-22)
* `-mv`, `-ms` で指定されていない原子が消えてしまうバグを修正した。

### Ver. 5.0 (2021-09-21)
* プログラム全体を書き直した。
	* 分子のトポロジーの処理を native の parmed に変更した。
	* 不要な独自モジュールを削除した。
	* Number モードが動作しないバグを修正した。
	* View mode と Extract mode をサブコマンド化した。

## Ver. 4.8 (2021-07-14)
* バージョン管理システムを mercurial から git へ変更した。

### Ver. 4.7 (2021-07-14)
* 公開した。
* PEP8 スタイルにした。
