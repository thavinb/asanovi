#!/usr/bin/env python

# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/asanovi samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,platform,fastq,fastq_1,fastq_2
    A1,nano,A1_LR.fastq.gz,A1_1.fastq.gz,A1_2.fastq.gz
    A2,nano,A2_LR.fastq.gz,,
    A3,,,A3_1.fastq.gz,A3_2.fastq.gz

    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 3
        HEADER = ["sample","platform","fastq","fastq_1", "fastq_2"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, platform, fastq, fastq_1, fastq_2 = lspl[: len(HEADER)]
            sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## DEBUG
            # debug = []
            # for i in [sample,platform,fastq,fastq_1,fastq_2]:
            #     if i:
            #         debug.append('TRUE')
            #     else:
            #         debug.append('FALSE')
            # print(sample,debug)
            # print(sample,[sample,platform,fastq,fastq_1,fastq_2])


            ## Check FastQ file extension
            for fq in [fastq, fastq_1, fastq_2]:
                if fq:
                    if fq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fq.endswith(".fastq.gz") and not fq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Check platform
            if platform not in ['nano','pacbio','']:
                print_error("Platform must specified either in 'nano','pacbio',or ''(blank)!")

            ## Auto-detect type of assembly
            sample_info = []  ## [ platform, fastq, fastq_1, fastq_2]
            if sample and not (platform and fastq) and fastq_1 and fastq_2:  ## Short reads assembly
                sample_info = ["shortread", platform, fastq, fastq_1, fastq_2]
                # print(sample,"shortread")
            elif sample and platform and fastq and not (fastq_1 and fastq_2):  ## Long read assembly
                sample_info = ["longread", platform, fastq, fastq_1, fastq_2]
                # print(sample,"longread")
            elif sample and platform and fastq and fastq_1 and fastq_2: ## Hybrid assembly
                sample_info = ["hybrid", platform, fastq, fastq_1, fastq_2]
                # print(sample,"hybrid")
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary
            ## { sample: [ platform, fastq, fastq_1, fastq_2 ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "method", "platform", "fastq", "fastq_1", "fastq_2"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error("Multiple runs of a sample must be of the same datatype!", "Sample: {}".format(sample))

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join(["{}".format(sample, idx + 1)] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
