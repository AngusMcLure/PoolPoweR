# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Upcoming]  

### Long-term
- [ ] Establish default values for functions [#11](https://github.com/AngusMcLure/PoolPoweR/issues/11)

### [v1.0.0]
- CRAN release
- Add pkgdown site (#20)  
- Implement codecov.io (#5)
- Add PoolPoweR-package.R (#18)  
- Configure and optimise GHA workflows

### [vX]
**Focus:** Enhancements  
- Refactor functions i.e. break down into smaller ones for readability (#26)

### [v0.1.3] - 2023-12-22
**Focus:** Documentation
- Revise documentation so majority @inheritParams fi_pool (#28)
- Rename interval (#13)  
- Update docs with cluster and non-cluster cases (#14) 
- Move descriptions to details  

### [v0.1.2] - 2023-12-15  
- Implement warnings for input checks (#3, #27)
- Include tests with expected inputs (#27, #34)  

### [v0.1.1] - 2023-12-08
- Add check_rho() to input_check() (#32)  

## [Available]  

### [v0.1.1] - 2023-12-07  

#### Added
- Class system `pool_strat_methods.R` for `catch_distribution()` (#15)  

#### Fixed
- Documentation for `_random()` functions and `pooling_strategies.R` (#31)
- Tests (#29, #30)
- Multiple local function definitions for `integrand` (#21)

### [v0.1.0] - 2023-12-05  

Period sampling functions.  

#### Added
- `fi_pool_cluster_random()` and `optimise_random_prevalence()` for evaluating
and designing surverys where the number of units caught at each site (or
overall) are random
- Machinery to make the above work (`catch_distributions.R`, `pooling_strategies.R`)

#### Fixed  
- Convergence of integrals with `fi_pool_cluster()` introduced by edge cases.  

#### Changed  
- `fi_pool_cluster()` integration will throw an error with extreme cases.  

### [v0.0.2] - 2023-12-05

Test coverage at 96.79%, but still requires considerable revision.

#### Added  
- Helper functions for arguments in `check_inputs.R`. Implemented in mainly
user-facing functions: `optimise_s_prevalence()`, `optimise_sN_prevalence()`,
`design_effect()` and `fi_pool_cluster()` (#2)
- Unit tests with reasonable examples from #11.  
- `test-coverage` and GHA (#5)

#### Changed  
- Error and warning handling for `fi_pool_cluster()` inputs

#### Removed
- Slower `fi_pool_cluster()` unit tests (cases where $s>1$ and $N>1$), these
also had extreme parameters.  

### [v0.0.1] MVP!  

#### Added  
- End-user-ready PoolPoweR package with core functions (#1, #4, #6, #7, #8,
	#10, #11)
- GitHub Action workflows `check-standard` and `test-coverage` implemented (#5)  

[upcoming]:
[0.0.2]:
[0.0.1]:
