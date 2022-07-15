# Changelog

All notable changes to **basinex** will be documented in this file.

## [0.2] - 2022-07

### Enhancements
- updated `report.out` with catchment area comparison and adjusted x, y of gauges ([#3](https://github.com/mhm-ufz/basinex/pull/3))
- `mask`, `gauge`, `gridfiles` and `ncfiles` are truly optional now ([#4](https://github.com/mhm-ufz/basinex/pull/4))

### Bugfixes
- `netcdf4` version needs to be `<1.6`, so we added a restriction to `setup.cfg` ([#3](https://github.com/mhm-ufz/basinex/pull/3))


## [0.1] - 2022-02

### Enhancements
- modern package structure: https://github.com/mhm-ufz/basinex
- added documentation: https://basinex.readthedocs.io
- made installable: `pip install basinex`
- added entry point for script usage
- added `latitude-size-correction` switch to yaml file for basin area correction with latlon coordinates
- added `--cwd` and `--version` to the basinex CLI

### Changes
- added ufz dependecies to src

### Bugfixes
- solved yaml warnings

[Unreleased]: https://github.com/mhm-ufz/basinex/compare/v0.1...HEAD
[0.2]: https://github.com/mhm-ufz/basinex/compare/v0.1...v0.2
[0.1]: https://github.com/mhm-ufz/basinex/releases/tag/0.1
