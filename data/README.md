# Data

All datasets are public. They are included here, unchanged except for the format noted, so that the scripts run without any download.

| File | Content | Source and licence |
|---|---|---|
| `barley/Durban2003Barley.csv` | Barley row–column trial: 272 genotypes × 2 replicates, 544 plots on a 16 × 34 grid (`row`, `bed`, `rep`, `gen`, `yield`) | Durban et al. (2003) *J. Agric. Biol. Environ. Stat.* 8, 48–66; distributed as `durban.rowcol` in the R package **agridat** (Wright, 2025; GPL) |
| `rice/Geno44k.txt` | Rice diversity panel: 413 accessions × 36,901 SNPs, coded 0/1/2 (NA = missing) | Zhao et al. (2011) *Nat. Commun.* 2, 467; [ricediversity.org](http://www.ricediversity.org/data/) (44K GWAS data set) |
| `rice/sativas413.map`, `rice/sativas413.fam` | PLINK map (chromosome, SNP, cM, bp; MSU6) and sample file for the same panel | as above |
| `rice/rice_phenotypes.csv` | 34 traits of the same 413 accessions (row names = `NSFTV_<id>`) | Zhao et al. (2011), as redistributed in `DT_rice` of the R package **enhancer** (GPL) |

Please cite the original publications when you reuse these data.

The PLINK `.ped` file is not needed by any script and is not tracked (see `.gitignore`).
