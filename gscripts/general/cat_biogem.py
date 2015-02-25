__author__ = 'gpratt'
import glob
import os
import argparse
from itertools import izip
import shutil


def warn_possibly_made(files):
    """
     Checks if all files are the same length, if they are don't do anything
     othewise warn and kill pasting, if ignore is set just warn, but don't kill
    """
    if len(set([len(fn) for fn in files])) != 1:
        print files
        raise ValueError("Files not all same length, possible malformation of folder, not writing to be safe")


def cat_strand(files, ext):
    if not files:
        raise ValueError("No files in directory")

    warn_possibly_made(files)
    #gets substring upto first difference in all selected files
    for x, file_chr in enumerate(izip(*files)):
        if len(set(file_chr)) != 1:
            break

    #get the file
    out_file = files[0][:x]

    #remove characters after last _
    out_file = "_".join(out_file.split("_")[:-1]) + ext

    #now that we have the outfile lets cat the actual files
    print "processing " + str(files) + " to " + out_file

    if not os.path.exists(out_file):
        with open(out_file, 'wb') as destination:
            for file_name in files:
                with open(file_name, 'rb') as file_obj:
                    shutil.copyfileobj(file_obj, destination)
    else:
        print "%s already made, not overwriting" % out_file


def cat_folder(dir_name, ext):
    r1 = glob.glob(os.path.join(dir_name, "*R1*%s" % ext))
    r2 = glob.glob(os.path.join(dir_name, "*R2*%s" % ext))
    cat_strand(r1, ext)
    cat_strand(r2, ext)


def walk_folders(directory, ext):
    for base, dir_names, files in os.walk(directory):
        try:
            cat_folder(base, ext)
        except ValueError as e:
            print e
            pass
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--recurse", action="store_true")
    parser.add_argument("-d", "--directory", dest="directory", default=".")
    parser.add_argument("-e", "--extension", dest="extension", default=".fastq.gz")
    parser.add_argument("")
    args = parser.parse_args()

    if not args.recurse:
        cat_folder(args.directory, args.extension)
    else:
        walk_folders(args.directory, args.extension)
