### 0.1.2
  * When sample is not specified, 'sample' is used as a placeholder (instead of an empty string).

### 0.2.0
* The report can now contain information about multiple samples.
  * `pairsqc.py` now takes sample name as an option (-s). The output files in the output report directory will have the sample name as file prefix.
  * `plot.r` also assumes sample names as prefix in the report directory (auto-detect).
  * Now d3 interactive plots are integrated with the Nozzle report. (https://github.com/SooLee/d3forNozzleR)

### 0.2.1
  * 'Sample' is used as a sample name place holder (instead of an empty string), when sample name is not specified.

### 0.2.2
  * `-M` option is introduced to work for non-human species (e.g. -M 8.2 for mouse)


