# This script runs unit tests to verify the filesystem initialized without errors

# Import modules
import os
from os.path import isfile, join
import numpy as np
import pandas as pd
from unit_tests import *
import subprocess


# Get all library ids
libraries = os.listdir("jobs")

# Unit test on info.txt
f = open("tests.out", "a")
lims = pd.read_csv("info.txt", sep = "\t", header = 0)
special_chars = "Do the sample ids contain special characters".format(special_char(lims))
f.writelines([special_chars])

# For each library, run the unit tests and save the results in tests.out
for lib in libraries:
    f.write("\nRunning unit tests on the library: {}\n".format(lib))
    
    # Unit tests on feature_reference.csv
    feat_ref = pd.read_csv("jobs/{}/feature_reference.csv".format(lib))
    f.write("\nTests on the file feature_reference.csv\n")
    ncols = "Does the file contain 6 columns? {}\n".format(feat_ref_ncols(feat_ref))
    header = "Is the header of the file correct? {}\n".format(feat_ref_header(feat_ref))
    columns = "Is the content of the columns the one expected? {}\n".format(feat_ref_columns(feat_ref))
    f.writelines([ncols, header, columns])

    # Unit tests on libraries.csv
    libs = pd.read_csv("jobs/{}/libraries.csv".format(lib))
    f.write("\nTests on the file libraries.csv\n")
    header_libs = "Is the header of the file correct? {}\n".format(libs_header(libs))
    nrows = "Does the file contain the correct number of rows? {}\n".format(libs_nrows(libs))
    lib_type = "Are the library types (cDNA/HTO) specified appropiately? {}\n".format(libs_library_type(libs))
    f.writelines([header_libs, nrows, lib_type])

    # Unit tests on fastq directories
    f.write("\nTests on fastq file directories\n")
    dirs = "Does the library have one fastq subdirectory per library type? {}\n".format(fastq_dirs(lib))
    files = "Does each fastq subdirectory has one fastq per read (R1/R2)? {}\n".format(fastq_files(lib))
    f.writelines([dirs, files])

# Unit tests on fastq file concatenation
f.write("\n\nTests on fastq file concatenation\n")
fastq = "Do the number of files of the original and concatenated fastq files match? {}\n".format(fastq_lines(lims))
f.write(fastq)
f.close()
