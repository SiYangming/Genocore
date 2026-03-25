# GenoCore : a simple and fast algorithm for core subset selection from large genotype datasets
Usage : Rscript run_genocore.R input_file -cv coverage -d difference -o result_file_name [--keep-all] [--seed SEED]

Genocore is available at https://github.com/lovemun/Genocore. Source code was written in R language and supported on windows and linux platform. 

## Requirement
- python : rst2pdf modules
  - install command in command line : pip install rst2pdf
- R : argparse package
  - install command in R : install.packages("argparse")

## Example

```shell
git clone https://github.com/lovemun/Genocore

cd Genocore

$ Rscript run_genocore.R wheat_subset.csv -cv 99 -d 0.001 -o example &
```
## --keep-all

默认情况下，程序会在满足以下任一条件时提前停止：

- Coverage >= -cv/--coverage
- Difference < -d/--delta

新增参数 --keep-all 后，不会因以上阈值提前停止，会继续按“贡献度/覆盖提升”的选择逻辑把所有样本完整排序输出（覆盖率达到 100% 后 Difference 通常为 0，但仍会继续选完剩余样本）。

### 用于比较开启/关闭 --keep-all 的示例

同一份输入分别跑两次（建议使用不同的 -o 前缀避免覆盖输出）：

$ Rscript run_genocore.R wheat_subset.csv -cv 99 -d 0.001 -o example_no_keep --seed 1

$ Rscript run_genocore.R wheat_subset.csv -cv 99 -d 0.001 -o example_keep --keep-all --seed 1

然后用 R 对两次结果做对比：

```r
no_keep_cov <- read.csv("example_no_keep_Coverage.csv", check.names = FALSE)
keep_cov <- read.csv("example_keep_Coverage.csv", check.names = FALSE)

cat("no_keep rows:", nrow(no_keep_cov), "\n")
cat("keep_all rows:", nrow(keep_cov), "\n")

stopifnot(nrow(keep_cov) >= nrow(no_keep_cov))

stopifnot(
  all(no_keep_cov$Sample_name == keep_cov$Sample_name[seq_len(nrow(no_keep_cov))]),
  all(abs(no_keep_cov$Coverage - keep_cov$Coverage[seq_len(nrow(no_keep_cov))]) < 1e-12)
)

no_keep_core <- read.csv("example_no_keep_Coreset.csv", check.names = FALSE)
keep_core <- read.csv("example_keep_Coreset.csv", check.names = FALSE)

stopifnot(ncol(keep_core) >= ncol(no_keep_core))
stopifnot(
  identical(names(no_keep_core), names(keep_core)[seq_len(ncol(no_keep_core))])
)

cat("Prefix一致性检查通过：不开启 --keep-all 的结果等于开启 --keep-all 结果的前缀。\n")
```

### 验证不开启 --keep-all 与之前结果是否一致

不开启 --keep-all 时，停止条件与历史版本保持一致（Coverage/Difference 触发即停止），因此在同一输入与参数下：

- example_no_keep_Coverage.csv 的行数与内容应与旧版本一致
- example_no_keep_Coreset.csv 的列数与列名顺序应与旧版本一致

## Contact

lovemun@kribb.re.kr

## Citation

Jeong S, Kim JY, Jeong SC, Kang ST, Moon JK, et al. (2017) GenoCore: A simple and fast algorithm for core subset selection from large genotype datasets. PLOS ONE 12(7): e0181420. https://doi.org/10.1371/journal.pone.0181420
