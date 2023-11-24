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
- [ ] What is the MVP required for a CRAN release?

### [v0.1.1]  
- [ ] Implement new functions  

### [v0.1.0]
**Focus:** Package
- [ ] Fix integrand issue  
- [ ] Implement codecov.io  
- [ ] Configure and optimise GHA workflows
- [ ] Refer to package and tidyverse guides for coherency and note rooms for improvement

### [v0.0.3]
**Focus:** Documentation
- [ ] Revise documentation so majority @inheritParams fi_pool
- [ ] Rename interval [#13](https://github.com/AngusMcLure/PoolPoweR/issues/13)
- [ ] Update docs with cluster and non-cluster cases
[#14](https://github.com/AngusMcLure/PoolPoweR/issues/14)
- [x] Rename real.scale to real_scale

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

## [Released]  

### [v0.0.1] MVP!  

#### Added  
- End-user-ready PoolPoweR package with core functions (#1, #4, #6, #7, #8,
	#10, #11)
- GitHub Action workflows `check-standard` and `test-coverage` implemented (#5)  

[unreleased]:
[0.0.1]:
