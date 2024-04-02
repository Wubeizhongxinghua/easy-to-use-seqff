# An easy-to-use version of seqff

> Seqff is tool for estimating fetal DNA fraction from https://obgyn.onlinelibrary.wiley.com/doi/full/10.1002/pd.4615

## How to use it?

+ You need to prepare:
	+ `.bam` file of a sample as input for estimating its fetal DNA fraction.
	+ All files listed in this repo., and put in the same directory.
	+ An R software with `Rsamtools`, `argparse` installed.

+ And then, just run the following code in terminal:

```shell
$ Rscript seqff.r -f input.bam -o output.txt
```
Where `input.bam` is file name of your `.bam` file, and `output.txt` is the file which will store the result of seqff.

## What have been revised?

+ You can input `.bam` files straightly after alignment instead of `.sam` files without header.
+ Sometimes, output of `enet` would be `NA` due to `NA` values in variable `bincounts`. We set them to zeros before calculating `enet` value.
+ Less arguments needed.
+ In some conditions, the chromosome names in `.bam` file do not contain string 'chr'. This script automatically detect it and add the 'chr' string (see [Issue #1](https://github.com/Wubeizhongxinghua/easy-to-use-seqff/issues/1)).
