# Developer resources and instructions  

**While developing**
- Regularly run `test()` or `test_active_file()`

**Before committing**
- Run `check()`, address any errors, warnings and notes

**Before pushing**
- Sync with origin/main
- Address any conflicts
- Run dependency graph
- Run `test-coverage`  

## Function dependency graph  

```r
#devtools::install_github("datastorm-open/DependenciesGraphs")
library(DependenciesGraphs)
load_all()
dep <- envirDependencies("package:PoolPoweR")
plot(dep)
```

![](https://raw.githubusercontent.com/AngusMcLure/PoolPoweR/main/dependencies/dep.png?token=GHSAT0AAAAAACICLPCIT5EJVC32Y7QQITYCZKJTYEA)  

## Test coverage  

```r
library(covr)
cov <- package_coverage()
write.csv(cov, "dev_resources/coverage_report.csv", quote = F, row.names = F)

# To show report markup in RStudio
report(cov)
```
