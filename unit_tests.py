# Contains one function per unit test.
# Checks if the output from 2-init.py is correct.
# Turn each bug into a unit test.

import numpy as np
import pandas as pd
import os
from os.path import isfile, join 


# Define unit tests as single functions
def feat_ref_ncols(df):
    """Tests if feature_reference.csv contains 6 columns.
    
    Args:
      df: pandas dataframe created from feature_reference.csv
    
    Returns:
      True if feature_reference.csv has 6 columns and False otherwise.
    """
    ncol = df.shape[1]
    if ncol == 6:
        return True
    else:
        return False

def feat_ref_header(df):
    """Tests if the header of feature_reference.csv is correct.
    
    Args:
      df: pandas dataframe created from feature_reference.csv. 
    
    Returns:
      True if feature_reference.csv has a correct header and False otherwise.
    """
    header = df.columns
    expected_header = ["id", "name", "read", "pattern", "sequence", "feature_type"]
    test = [i == j for i, j in zip(header, expected_header)]
    if test.count(True) == 6:
       return True
    else:
       return False

def feat_ref_columns(df):
    """Tests if the content of the columns of feature_reference.csv is correct.
    
    Args:
      df: pandas dataframe created from feature_reference.csv.
    
    Returns:
      True if the content of feature_reference.csv is correct, False otherwise.
    """
    # Are the first two columns equal?
    col1 = np.array(df[["id"]])
    col2 = np.array(df[["name"]])
    if np.sum(col1 == col2) != df.shape[0]:
        return False
    
    # Is the hashtag located in read 2?
    col3 = np.array(df[["read"]])
    if np.sum(col3 == "R2") != df.shape[0]:
        return False
    
    # Does the 4rt column contain the pattern ^(BC)?
    col4 = np.array(df[["pattern"]])
    if np.sum(col4 == "^(BC)") != df.shape[0]:
        return False

    # Does the last column contain "Antibody Capture"?
    col6 = np.array(df[["feature_type"]])
    if np.sum(col6 == "Antibody Capture") != df.shape[0]:
        return False
    else:
        return True

def libs_header(df):
    """Tests if the header of libraries.csv is correct.
    
    Args:
      df: pandas dataframe created from libraries.csv.
    
    Returns:
      True if libraries.csv has a correct header, False otherwise.
    """
    if df.columns[0] == "fastqs" and df.columns[1] == "sample" and df.columns[2] == "library_type":
        return True
    else:
        return False

def libs_nrows(df):
    """Tests if libraries.csv has 2 rows (without counting header).
    
    Args:
      df: pandas dataframe created from libraries.csv.
    
    Returns:
      True if libraries.csv has 2 rows, False otherwise.
    """
    if df.shape[0] == 2:
        return True
    else:
        return False

def libs_library_type(df):
    """Tests if the correspondence between library type (cDNA/HTO) and fastqs path is correct.
    
    Args:
      df: pandas dataframe created from libraries.csv.
    
    Returns:
      True if the correspondence is correct, False otherwise.
    """
    fastqs = df["fastqs"].values
    lib_types = df["library_type"].values
    lib_type_dict = {"HTO":"Antibody Capture", "cDNA":"Gene Expression"}
    fastq_dirs = [fastqs[x].split("/")[-1] for x in [0,1]]
    if lib_type_dict[fastq_dirs[0]] == lib_types[0] and lib_type_dict[fastq_dirs[1]] == lib_types[1]:
        return True
    else:
        return False

def fastq_dirs(library):
    """Tests if a library has two separate fastq subdirectories for each library type (cDNA and HTO) 
    
    Args:
      library: string specifying the library id.
    
    Returns:
      True if the fastq subdirectories exist and are in the right path, false otherwise.
    """
    fastq_dirs = os.listdir("jobs/{}/fastq/".format(library))
    if len(fastq_dirs) == 2 and "cDNA" in fastq_dirs and "HTO" in fastq_dirs:
        return True
    else:
        return False

def fastq_files(library):
    """Tests if there are two fastq files (read 1 and read2) for each library type (cDNA and HTO
    
    Args:
      library: string specifying the library id.
    
    Returns:
      True if there are two fastq files per library, false otherwise.
    """
    for lib_type in ["cDNA", "HTO"]:
        path = "jobs/{}/fastq/{}".format(library, lib_type)
        files = [f for f in os.listdir(path) if isfile(join(path, f))]
        if len(files) != 2:
            return False
        sufix_r1 = "_S1_L001_R1_001.fastq.gz"
        sufix_r2 = "_S1_L001_R2_001.fastq.gz"
        if not(sufix_r1 in files[0] or sufix_r1 in files[1]):
            return False
        if not(sufix_r2 in files[0] or sufix_r2 in files[1]):
            return False
    return True

def fastq_lines(lims):
    """Tests if the number of lines of the concatenated fastqs is the expected.
    
    Args:
      lims: pandas dataframe with the information of the lims for that subproject.
    
    Returns:
      True if the fastqs were concatenated properly, false otherwise.
    """
    import subprocess
    for lib_type in ["cDNA", "HTO"]:
        if lib_type == "HTO":
            lims_sub = lims[["Illumina" in x for x in lims["index"]]]
        else:
            lims_sub = lims[["Illumina" not in x for x in lims["index"]]]
        
        # Search fastqs in /project from scratch and compute the total number of lines
        fc_lane_ind = [(lims_sub["flowcell"][x], lims_sub["lane"][x], lims_sub["index"][x]) for x in lims_sub.index]
        fastq_project = "/project/production/fastq"
        fastq_original = ["{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(fastq_project, fc, lane, fc, lane, index) for fc, lane, index in fc_lane_ind]
        nlines_list = []
        for path in fastq_original:
            cmd = "zcat {} | wc -l".format(path)
            ps = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            output = ps.communicate()
            nlines = int(output[0].decode("utf-8").strip("\n"))
            nlines_list.append(nlines)
        nlines_original = sum(nlines_list)
        
        # Calculate total number of lines concatenated fastq files
        fastq_concatenated = ["jobs/{}/fastq/{}/".format(x, lib_type) for x in os.listdir("jobs")]
        fastq_concatenated = ["{}/{}".format(x, os.listdir(x)[0]) for x in fastq_concatenated]
        nlines_list = []
        for path in fastq_concatenated:
            cmd = "zcat {} | wc -l".format(path)
            ps = subprocess.Popen(cmd, shell = True, stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
            output = ps.communicate()
            nlines = int(output[0].decode("utf-8").strip("\n"))
            nlines_list.append(nlines)
        nlines_concatenated = sum(nlines_list)
        
        # Compare and return
        if nlines_original == nlines_concatenated:
            return True
        else:
            return False

def fastq_exist(lims):
    """Tests if all the fastqs exist.
    
    Args:
      lims: pandas dataframe with the information of the lims for that subproject.
    
    Returns:
      True if the fastqs exist, false otherwise.
    """
    fastq_path = "/project/production/fastq"
    fc_lane_ind_list = [(lims["flowcell"][x], lims["lane"][x], lims["index"][x]) for x in lims.index]
    fastq_path_list_r1 = ["{}/{}/{}/fastq/{}_{}_{}_1.fastq.gz".format(fastq_path, fc, lane, fc, lane, ind) for fc, lane, ind in fc_lane_ind_list]
    fastq_path_list_r2 = ["{}/{}/{}/fastq/{}_{}_{}_2.fastq.gz".format(fastq_path, fc, lane, fc, lane, ind) for fc, lane, ind in fc_lane_ind_list]
    test_r1 = np.sum([os.path.isfile(x) for x in fastq_path_list_r1]) == lims.shape[0]
    test_r2 = np.sum([os.path.isfile(x) for x in fastq_path_list_r2]) == lims.shape[0]
    if test_r1 and test_r2:
        return True
    else:
        return False

def special_char(lims):
    """Tests if sample ids do not contain special characters
    
    Args:
      lims: pandas dataframe with the information of the lims for that subproject.
    
    Returns:
      True if the sample ids do not contain special characters.
    """
    chars = string.punctuation.replace("-", "").replace("_", " ")
    tests = [any(char in chars for char in x) for x in lims["SampleName"]]
    if np.sum(tests) == 0:
        return True
    else:
        return False
