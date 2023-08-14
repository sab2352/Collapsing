#!/usr/bin/env python

"""! @brief Python program to review samples removed to create ped.txt and sample_file, and identify any case samples removed"""

import os
import argparse
from datetime import datetime

def main():
    """! Main program entry."""
    parser = argparse.ArgumentParser(description='Program to review samples removed to create ped.txt and sample_file')
    parser.add_argument('--caselist', type=str, help='file path for case list', required=True)
    parser.add_argument('--ped_txt', type=str, help='file path for ped.txt', required=True)
    parser.add_argument('--sample_file_check', type=str, help='file path for sample_file', required=True)
    parser.add_argument('--out', type=str, help='output directory', required=True)
    args = parser.parse_args()
    
    #check if output dir exists, if not make one
    pathExist = os.path.exists(args.out)
    if not pathExist:
        os.makedirs(args.out)
        
    #open needed files
    ped_txt = open(args.ped_txt)
    cases = open(args.caselist)
    sample_file = open(args.sample_file_check)
    
    #open/create and add time/date to all file outputs
    now = datetime.now()
    case_missing_diff_ped_existing = open("{0}/{1}_case_samples_removed_diff_ped_existing.txt".format(args.out, now.strftime("%y-%m-%d_%H-%M-%S")), "w")
    case_missing_ped = open("{0}/{1}_case_samples_removed_ped.txt".format(args.out, now.strftime("%y-%m-%d_%H-%M-%S")), "w")
    diff_ped_existing = open("{0}/{1}_diff_ped_existing.txt".format(args.out, now.strftime("%y-%m-%d_%H-%M-%S")), "w")
    log = open("{0}/{1}_sample_check_log.txt".format(args.out, now.strftime("%y-%m-%d_%H-%M-%S")), "a")  
    
    ped_list = make_list(ped_txt, True)
    existing_list = make_list(sample_file, True)
    case_list = make_list(cases, False)

    ped_txt.close()
    cases.close()
    sample_file.close()

    #compare files and see what samples have been removed
    diff_list = compare_ped_existing(ped_list, existing_list, diff_ped_existing, log)
    compare_ped_cases(ped_list, case_list, case_missing_ped, log)
    compare_diff_cases(diff_list, case_list, case_missing_diff_ped_existing, log)
    
    log.close()

#make input files into lists
def make_list(file_name, atav_out):
    """! Generate a list from the input file

    @param file_name File: Input file to create list
    @param atav_out Boolean: State whether file is the output of atav

    @return The generated list 
    """
    sample_list = []
    if atav_out:
        for line in file_name:
            sampleID = line.split()
            sample_list.append(sampleID[0].casefold())
    else:
        for line in file_name:
            sample_list.append(line.strip().casefold())
    
    return sample_list

#compare and see what was removed frome ped.txt to create sample_file    
def compare_ped_existing(ped_list, existing_list, diff_ped_existing, log):
    """! Compare ped.txt to sample_existing.txt and identify what samples were removed 

    @param ped_list List: List of samples in ped.txt   
    @param existing_list List: List of samples in existing.txt
    @param diff_ped_existing File: Output file for sample difference between lists
    @param log File: Output log file

    @return: List of samples what were removed from ped.txt to create sample_file
    """
    num_missing = 0
    diff_list = []
    
    log.write("{} \n".format(datetime.now()))
    log.write("Now reviewing changes made to ped.txt to create sample_file_existing... \n")
    log.write("Comparing samples present in ped.txt with samples present in existing.sample.txt \n")
    log.write("Results in diff_ped_existing.txt\n")
    print("Now reviewing changes made to ped.txt to create sample_file_existing...")
    
    for sample in ped_list:
        if sample not in existing_list:
            diff_ped_existing.write(sample + "\n")
            num_missing += 1
            diff_list.append(sample)
    diff_ped_existing.close()
    
    print("Samples removed from ped.txt:", num_missing)
    log.write("{} \n".format(datetime.now()))
    log.write("Number of samples removed from ped.txt is {} \n".format(num_missing))

    return diff_list

#compare and see what cases samples were removed from ped.txt
def compare_ped_cases(ped_list, case_list, case_missing_ped, log):
    """! Compare ped.txt to case list and identify what case samples were removed 

    @param ped_list List: List of samples in ped.txt   
    @param case_list List: List of case samples
    @param case_missing_ped File: Output file for sample difference between lists
    @param log File: Output log file  
    """
    num_missing = 0
    
    log.write("{} \n".format(datetime.now()))
    log.write("Comparing case sample list to ped.txt and identifying case samples removed to create ped.txt... \n")
    log.write("Results in case_samples_removed_ped.txt\n")
    print("Identifying case samples removed to create ped.txt")
    case_missing_ped.write("List of case samples not in ped.txt: \n")
      
    for sample in case_list:
        if sample not in ped_list:
            case_missing_ped.write(sample + "\n")
            num_missing += 1
    case_missing_ped.close()
    
    print("Number of case samples not in ped.txt is", num_missing)
    log.write("{} \n".format(datetime.now()))
    log.write("Number of case samples not in ped.txt is {} \n".format(num_missing))
    
#compare and see which case samples were among ones removed to create sample_file  
def compare_diff_cases(diff_list, case_list, case_missing_existing, log):
    """! Compare case list to difference between ped.txt and coverage_exisiting and identify what case samples were removed 

    @param diff_list List: List of samples removed from ped.txt to create sample_file  
    @param case_list List: List of case samples
    @param case_missing_existing File: Output file for sample difference between lists
    @param log File: Output log file   
    """
    num_missing = 0
    
    log.write("{} \n".format(datetime.now()))
    log.write("Identifying case samples removed from ped.txt to create existing.sample.txt... \n")
    log.write("Comparing case samples present in ped.txt with case samples present in existing.sample.txt \n")
    log.write("Results in case_samples_removed_diff_ped_existing\n")
    print("Identifying case samples removed from ped.txt to create sample_file")
    case_missing_existing.write("List of case samples removed from ped.txt to create sample_file: \n")
    
    for sample in case_list:
        if sample in diff_list:
            case_missing_existing.write(sample + "\n")
            num_missing += 1
    case_missing_existing.close()
    
    print("Number of case samples among ones removed from ped.txt is", num_missing)
    log.write("{} \n".format(datetime.now()))
    log.write("Number of case samples among ones removed from ped.txt to create sample_file is {} \n".format(num_missing))
    

if __name__=="__main__":
    main()
