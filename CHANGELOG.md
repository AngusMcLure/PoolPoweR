# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Upcoming]  

### Long-term
- [ ] Establish default values for functions [#11](https://github.com/AngusMcLure/PoolPoweR/issues/11)
- [ ] Create full class system for catch_distribution [#15](https://github.com/AngusMcLure/PoolPoweR/issues/15)

### [v1.0.0]
- CRAN release
- Add pkgdown site #20  
- [ ] What is the MVP required for a CRAN release?

### [v0.1.1]  
- [ ] Implement new functions  

### [v0.1.0]
**Focus:** Package
- [ ] Fix integrand issue  
- [ ] Implement codecov.io  
- [ ] Refactor functions (i.e. break down into smaller ones for readability)
- [ ] Configure and optimise GHA workflows

### [v0.0.3]
**Focus:** Documentation
- [ ] Revise documentation so majority @inheritParams fi_pool
- [ ] Rename interval [#13](https://github.com/AngusMcLure/PoolPoweR/issues/13)
- [ ] Update docs with cluster and non-cluster cases
[#14](https://github.com/AngusMcLure/PoolPoweR/issues/14)
- [ ] Add PoolPoweR-package.R (#18)  
- [ ] Move descriptions to details  
- [x] Rename real.scale to real_scale

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

## [Released]  

### [v0.0.1] MVP!  

#### Added  
- End-user-ready PoolPoweR package with core functions (#1, #4, #6, #7, #8,
	#10, #11)
- GitHub Action workflows `check-standard` and `test-coverage` implemented (#5)  

[unreleased]:
[0.0.1]:
