# An easy-to-use version of seqff

> Seqff is tool for estimating fetal DNA fraction from https://obgyn.onlinelibrary.wiley.com/doi/full/10.1002/pd.4615

## How to use it?

+ You need to prepare:
	+ `.bam` file of a sample as input for estimating its fetal DNA fraction.
	+ All files listed in this repo., and put in the same directory.
	+ An R software with `Rsamtools`, `argparse`, `dplyr`, `stringr` installed.

+ And then, just run the following code in the terminal:

```shell
$ Rscript seqff.r -f input.bam -o output.txt
```
Where `input.bam` is file name of your `.bam` file, and `output.txt` is the file which will store the result of seqff.

If your `.bam` file was aligned to hg38 genome instead of the default hg19 genome, you can run by:
**When using hg38 mode, you shall read the disclaimer (on the end of the page) first!**
```shell
$ Rscript seqff.r -f input.bam -o output.txt --hg38
```

## What have been revised?

+ You can input `.bam` files straightly after alignment instead of `.sam` files without header.
+ Sometimes, output of `enet` would be `NA` due to `NA` values in variable `bincounts`. We set them to zeros before calculating `enet` value.
+ Less arguments needed.
+ In some conditions, the chromosome names in `.bam` file do not contain string 'chr'. This script automatically detects them and adds the 'chr' strings (see [Issue #1](https://github.com/Wubeizhongxinghua/easy-to-use-seqff/issues/1)).
+ [Hg38 mode](https://github.com/Wubeizhongxinghua/easy-to-use-seqff/issues/2) have been implemented by setting `--hg38` (see help document for details). This mode is achieved by [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver). 
> Disclaimer: When using `liftOver` for conversion, 8.29% of bins cannot correspond to hg38. Therefore, the author cannot guarantee that the hg38 mode necessarily reflects the ideal performance of the model. The author tested a small number of samples (aligned to hg38 genome, using hg38 and hg19 modes respectively) and found that the predicted results remained almost unchanged.
