from pyflow import WorkflowRunner
from optparse import OptionParser
import glob
import os
import config
import cPickle
import random
import utils

# this workflow checks the hudson alpha fastqs

def get_params():
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-i","--input_dir",
                      dest="source",
                      default=None,
                      help="input directory")
    parser.add_option("-j","--jobs",
                      dest="jobs",
                      default=10,
                      help="maximum number of jobs to submit concurrently")
    parser.add_option("-o","--output_dir",
                      dest="dest",
                      default=None,
                      help="output directory")
    parser.add_option("--slids",
                      dest="slids",
                      default=None,
                      help="only analyze specified samples by SLID, e,g. SL32616.  Use a comma to analyze multiple SLIDs.")

    parser.add_option("-p","--plate",
                      dest="plate",
                      default=None,
                      help="if specified, only analyze samples from plate.")
    parser.add_option("-l","--local",
                      dest="local",
                      action='store_true',
                      default=False,
                      help="run on cluster")

    (options,args) = parser.parse_args()
    if None in [options.source,options.dest]:
        print "error.  must specify source and destination directories."
        parser.print_help()
        sys.exit(-1)
    options.jobs = int(options.jobs)
    if options.slids != None:
        options.slids = set(options.slids.split(","))

    return options,args


class PyflowMD5Checker(WorkflowRunner):
    def __init__(self, input_directory, output_directory,plate,slids):
        self.input_directory = input_directory
        self.output_directory = output_directory
        self.plate = plate
        self.slids = slids
        print "getting slids"
        self.slid_mapper = utils.get_slids_from_dir(self.input_directory,slids=self.slids)
        print self.slid_mapper.keys()
        print self.slid_mapper[self.slid_mapper.keys()[0]].keys()

        
    def workflow(self):
        tasknum = 0
        plates = self.slid_mapper.keys()
        if not self.plate == None:
            plates = [self.plate]

        for plate in plates:
            for slid in self.slid_mapper[plate]:
                input_dir = os.path.join(self.input_directory,plate)
                output_dir = os.path.join(self.output_directory,plate)
                taskname = "md5_%i" % (tasknum)
                tasknum += 1
                dependencies=()
                python = config.anaconda_python
                script = os.path.join(config.python_src_dir,
                                      "single_slid_md5sum_check.py")
            
                command = "%s %s -s %s -i %s -o %s" % (python,
                                                       script,
                                                       slid,
                                                       input_dir,
                                                       output_dir)
                self.addTask(taskname,
                             command=command,
                             dependencies=dependencies)



if __name__=='__main__':
    options, args = get_params()
    if not os.path.exists(options.dest):
        os.makedirs(options.dest)

    fs = PyflowMD5Checker(options.source,
                          options.dest,
                          options.plate,
                          options.slids)

    mode = None
    if options.local:
        mode = "local"
    else:
        mode = "torque"

    print "executing md5sum checker workflow. mode:",mode
    fs.run(mode=mode,
           dataDirRoot=options.dest,
           nCores=options.jobs)

           
           
        




