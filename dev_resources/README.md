# Function dependencies  

Code to generate the function dependency graph:  
```r
devtools::install_github("datastorm-open/DependenciesGraphs")
library(DependenciesGraphs)
load_all()
dep <- envirDependencies("package:PoolPoweR")
plot(dep)
```

![](https://raw.githubusercontent.com/AngusMcLure/PoolPoweR/main/dependencies/dep.png?token=GHSAT0AAAAAACICLPCIT5EJVC32Y7QQITYCZKJTYEA)
