# tcga-vcf-reheader-tool

Tool to read a TCGA Variant Call Format (VCF) file and output an equivalent file with a different header

## Requirements

* Python 2.7
* PyYAML

## Usage

    tcga-vcf-reheader.py input_file_path output_file_path parameter_file_path

* input_file_path: the VCF to read
* output_file_path: the VCF to write
* parameter_file_path: a YAML file describing the additions/changes to the header

Currently only changes to the `##SAMPLE` lines are supported. Any existing `##SAMPLE` lines will be replaced by the ones specified in the YAML file.

**WARNING:** This tool will **completely replace** any existing ##SAMPLE lines. If the parameter file names different samples than were used in the VCF, the resulting output will be a misleading chimeric monster.

## Contributors

* Kyle Ellrott <kellrott@soe.ucsc.edu>
* Walker Hale <hale@bcm.edu>
* Dave Larson <dlarson@genome.wustl.edu>
* Lee Lichtenstein <lichtens@broadinstitute.org>
# tcgavcf-tool
