# Function dependencies  

Code to generate the function dependency graph:  
```r
devtools::install_github("datastorm-open/DependenciesGraphs")
library(DependenciesGraphs)
load_all()
dep <- envirDependencies("package:PoolPoweR")
plot(dep)
```

![](https://github.com/AngusMcLure/PoolPoweR/dependencies/dep.png)
