from pyflow import WorkflowRunner
from optparse import OptionParser
import glob
import os
# this is a simple test script that I am using to get pyflow working on TSCC

def get_params():
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage)
    parser.add_option("-s","--source_dir",
                      dest="source",
                      default=None,
                      help="input directory")
    parser.add_option("-d","--dest_dir",
                      dest="dest",
                      default=None,
                      help="output directory")
    parser.add_option("-c","--cluster",
                      dest="local",
                      action='store_false',
                      default=True,
                      help="run on cluster")

    (options,args) = parser.parse_args()
    if None in [options.source,options.dest]:
        print "error.  must specify source and destination directories."
        parser.print_help()
        sys.exit(-1)
    return options,args


# just does "ls -l" on each file in a directory; this simple workflow is
# intended to get pyflow working on tscc
class FileSizer(WorkflowRunner):
    def __init__(self, input_directory, output_directory):
        self.input_directory = input_directory
        self.output_directory = output_directory

    def workflow(self):
        files = glob.glob(os.path.join(self.input_directory,"*"))
        mid_file = files[len(files)/2]
        task_idx = 1
        for result_filename in files:
            filename = os.path.basename(result_filename)
            task_name = "test_task_%s" % (task_idx)
            introduce_dependency = False
            # the second half can't go until the first half
            # are done
            if task_idx == len(files)/2:
                dependency_name = task_name
            if task_idx > len(files)/2:
                introduce_dependency = True
            dependencies = ()
            if introduce_dependency:
                dependencies = ( dependency_name )
            input_name = os.path.join(self.input_directory,filename)
            output_name = os.path.join(self.output_directory,filename+".lsfilesize")
            command_line = "/bin/ls -lh %s > %s" % (input_name,
                                                output_name)
            self.addTask(task_name,
                         command=command_line,
                         dependencies=dependencies)
                         
            task_idx = task_idx + 1            
            
                
    
if __name__=='__main__':
    options, args = get_params()
    if not os.path.exists(options.dest):
        os.makedirs(options.dest)
    fs = FileSizer(options.source,
                   options.dest)
    if options.local:
        print "local."
        fs.run(mode="local")
    else:
        print "not local!"
        fs.run(mode="torque")
           
           
        




