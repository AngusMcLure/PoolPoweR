# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### [v0.0.1] - 2023-11-24
**Focus:** Package
- Minimum viable package online 
- README.md updated i.e. installation 
- GitHub actions implemented `check-standard`
- No errors or warnings in R CMD Check and GHA workflows
- [ ] Fix undocumented warning (cost_fi, cost_fi_unclustered, all utils.R)
	- [ ] Move utils functions (clogclog, cloglog_inv, mu_sigma_linknorm) 
	to fi_pool_cluster as only called once each. Will depend on other
	functions to be added
	- [ ] Move cost_fi and cost_fi_cluster to optimise_s_prevalence
	- [ ] or, simply document the helper functions

### [v0.0.2] - 2023-12-01
**Focus:** Tests
- [ ] Add helper function for input variable checks
[#2](https://github.com/AngusMcLure/PoolPoweR/issues/2)
- [ ] Throw warnings for some input variable ranges
[#3](https://github.com/AngusMcLure/PoolPoweR/issues/3)
	- Note: these parts currently lack test coverage in user-level functions
- [ ] Update unit tests with examples added in
[#11](https://github.com/AngusMcLure/PoolPoweR/issues/11)
- [ ] Add GHA for test coverage
[#5](https://github.com/AngusMcLure/PoolPoweR/issues/5)
- [ ] Add tests for uncovered lines
- [ ] Speed up tests (some slow ones fo fisher_info.R)

### [v0.0.3]
**Focus:** Documentation
- [ ] Revise documentation so majority @inheritParams fi_pool
- [ ] Rename interval [#13](https://github.com/AngusMcLure/PoolPoweR/issues/13)
- [ ] Update docs with cluster and non-cluster cases
[#14](https://github.com/AngusMcLure/PoolPoweR/issues/14)
- [x] Rename real.scale to real_scale

### [v0.0.4]
**Focus:** Package
- [ ] Refer to package and tidyverse guides for coherency and note rooms for improvement

### [v1.0.0]
- CRAN release
- [ ] What is the MVP for a CRAN release?

### Long-term
- [ ] Establish default values for functions [#11](https://github.com/AngusMcLure/PoolPoweR/issues/11)
