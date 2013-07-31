
from collections import defaultdict
import re
import math
import subprocess
from subprocess import PIPE

class Submitter:
  """
  Class that will customize and submit shell scripts
  to the job scheduler
  """

  def __init__(self, **kwargs):
    """
    Constructor method, will initialize class attributes to
    passed keyword arguements and values.
    """
    self.data = kwargs

  def set_value(self, **kwargs):
    """ 
    Set class attributes by passing keyword arguement and values
    """
    for key in kwargs:
      self.data[key] = kwargs[key]

  def add_wait(self, wait_ID):
    """ 
    Add passed job ID to list of jobs for this job submission to
    wait for. Can be called multiple times.
    """
    if 'wait_for' not in  self.data:
      self.data['wait_for'] = []
    
    self.data['wait_for'].append(str(wait_ID))

  def add_resource(self, kw, value):
    """
    Add passed keyword and value to a list of attributes that
    will be passed to the scheduler
    """
    if 'addtl_resc' not in self.data:
      self.data['addtl_resc'] = defaultdict(list)

    self.data['addtl_resc'][kw].append(value)

  def write_sh(self, **kwargs):
    """
    Create a file that is the shell script to be submitted to the
    job scheduler. 

    keyword arguement attributes that can be passed:
    array: Oolite only. distribute jobs to array of jobs if True. 
    chunks: if array is True, Int number of commands per job
    submit: submit to scheduler if True, only write SH file if False or None
    queue: TSCC only. name of queue to submit jobs to. defaults to 'home'
    walltime: TSCC only. specify a walltime for this job as HH:MM:SS string. defaults to 72h
    nodes: TSCC only. specify number of nodes to use. defaults to 1
    ppn: TSCC only. specify number of processors per node. defaults to 1

    method returns job ID assigned by scheduler for this job if submit is True,
    returns 0 if only writing SH file.

    Requires these class attributes:
    queue_type: Scheduler type, either SGE or PBS
    sh_file: name of shell file to write
    command_list: list of commands to be performed for this job
    job_name: name of this job

    Optional class attributes:
    wait_for: list of ID's that this job will wait for before starting
    addtl_resc: keyword-value pairs of scheduler attributes, set with add_resource()
    """
    required_keys = "queue_type sh_file command_list job_name".split()
    use_array = False
    chunks = 1
    submit = False
    ret_val = 0

    if 'array' in kwargs:
      use_array = kwargs['array']

    if 'chunks' in kwargs:
      if kwargs['chunks']:
        chunks = kwargs['chunks']
    
    if 'submit' in kwargs:
      submit = kwargs['submit']

    if 'walltime' in kwargs:
      walltime = kwargs['walltime']
    else:
      walltime = "72:00:00"

    if 'nodes' in kwargs:
      nodes = kwargs['nodes']
    else:
      nodes = 1

    if 'ppn' in kwargs:
      ppn = kwargs['ppn']
    else:
      ppn = 1

    if 'account' in kwargs:
      account = kwargs['account']
    else:
      account = 'yeo-group'

    if 'queue' in kwargs:
      queue = kwargs['queue']
    else:
      queue = 'home'

    for keys in required_keys:
      if not keys in self.data:
        print "missing required key: "+str(keys)
        return

      if not self.data[keys]:
        print "missing value for required key: "+str(keys)
      	return

    
    sf = open(self.data['sh_file'], 'w')

    if (self.data['queue_type'] == 'SGE'):

      sf.write("#!/bin/sh\n")
      sf.write("#$ -N "+self.data['job_name']+"\n")
      sf.write("#$ -o "+self.data['job_name']+".out\n")
      sf.write("#$ -e "+self.data['job_name']+".err\n")
      sf.write("#$ -V\n")
      sf.write("#$ -S /bin/sh\n")
      sf.write("#$ -cwd\n")

      if 'wait_for' in self.data:
        if self.data['wait_for']:	  
          sf.write("#$ -hold_jid "+" ".join(self.data['wait_for'])+"\n")

      if 'addtl_resc' in self.data:
        if self.data['addtl_resc']:
          for key in self.data['addtl_resc']:
            for value in self.data['addtl_resc'][key]:
              sf.write("#$ "+key+" "+value+"\n")

      if use_array:
        print "use array"
        number_jobs = math.ceil(len(self.data['command_list'])/int(chunks))

      else:
        for command in self.data['command_list']:
          sf.write(str(command) + "\n")

      sf.close()
      if submit:

        p = subprocess.Popen(["qsub", self.data['sh_file']], stdout=PIPE)
      	output = p.communicate()[0].strip()
        ret_val = re.findall(r'\d+', output)[0]

      return ret_val

    elif (self.data['queue_type'] == 'PBS'):

      sf.write("#!/bin/sh\n")
      sf.write("#PBS -N "+self.data['job_name']+"\n")
      sf.write("#PBS -o "+self.data['job_name']+".out\n")
      sf.write("#PBS -e "+self.data['job_name']+".err\n")
      sf.write("#PBS -V\n")
      sf.write("#PBS -l walltime={}\n".format(walltime))
      sf.write("#PBS -l nodes={}:ppn={}\n".\
        format(str(nodes), str(ppn)))
      sf.write("#PBS -A {}\n".format(account))
      sf.write("#PBS -q {}\n".format(queue))
      sf.write("cd $PBS_O_WORKDIR\n")

      if 'wait_for' in self.data:
        if self.data['wait_for']:   
          sf.write("#$ -hold_jid "+" ".join(self.data['wait_for'])+"\n")

      if 'addtl_resc' in self.data:
        if self.data['addtl_resc']:
          for key in self.data['addtl_resc']:
            for value in self.data['addtl_resc'][key]:
              sf.write("#$ "+key+" "+value+"\n")

      if use_array:
        print "use array"
        number_jobs = math.ceil(len(self.data['command_list'])/int(chunks))

      else:
        for command in self.data['command_list']:
          sf.write(str(command) + "\n")

      sf.close()
      if submit:

        if 'wait_for' in self.data:
          wait_cmd = '-W depend=\"afterok:\"'+':'.join(self.data['wait_for'])
          p = subprocess.Popen(["qsub", wait_cmd, self.data['sh_file']], stdout=PIPE)

        else:
          p = subprocess.Popen(["qsub", self.data['sh_file']], stdout=PIPE)

        output = p.communicate()[0].strip()
        ret_val = re.findall(r'\d+', output)[0]

      return ret_val

# to do: separate general queue methods from queue_type specifc methods.
# wrap each queue-type into a separate sub-class that inherits methods from Submitter


acceptableQueueTypes = set(['SGE', "PBS"])


class InteractiveSubmitter(object): #use ClassName(object) to use new python class
  """
  Class that will customize and submit shell scripts
  to the job scheduler

  #full invocation
  qS = Submitter()
  qS['queueType'] = 'SGE'
  qS['jobName'] = 'testJob'
  qS.add_Q_resource('-l', 'h_vmem=60G')
  qS.add_Q_resource('-l', 'bigmem')
  qS.add_command('echo started')
  qS.add_command('sleep 5')  
  qS.add_command('echo finished')
  qS.write_sh(shFile='testShFile.sh')
  qS.submit()

  #lazy invocation
  qS = Submitter()
  qS.submit(['echo \"this started is done\"', 'sleep 5', 'echo \"this job is done\"'],
  job_name='testJob', shFile='test.sh')
  
  #by default oolite jobs are submitted as arrays and run in parallel
  #maximum concurrent array job limit is auto-set to 15, modify this with:
  qS.maxJobs = 10
  
  #set useArray=False when calling submit() or write_sh() to avoid this.
  
  """
  
  def __init__(self, printHost=True, **kwargs):
    """
    
    Constructor method, will initialize class attributes to
    passed keyword arguements and values.

    by default printHost == True and the first line of stdout
    will be the name of the machine where a job is running.
    
    """
    self.maxJobs= None
    self.data = kwargs    
    self.data['shWasWritten']=False #has an sh script been created yet?

    if 'addtlResc' not in self.data:
      self.data['addtlResc'] = defaultdict(list)
    if 'waitFor' not in  self.data:
      self.data['waitFor'] = set()

    self.data['useArray'] = False
    self.data['commandList'] = list()
    self.data['sentJobs'] = list() #jobs that have been submitted
    try: #if you know the environment, get information for lazy people like me.
      self.auto_detect_queue()
      print "auto-detected queue settings"
    except:
      pass


  def auto_detect_queue(self):
    """ set parameters smart-ly"""

    hostname = subprocess.Popen('hostname', stdout=subprocess.PIPE).communicate()[0].strip()

    if ("oolite" in hostname) or ("compute" in hostname):
      self.data['queueType'] = "SGE"
      self.data['useArray'] = True
      self.add_Q_resource("-l", 'bigmem')
      self.add_Q_resource("-l", 'h_vmem=16G')
      self.maxJobs=15
    elif ("tscc" in hostname):
      self.data['queueType'] = "PBS"
      self.data['useArray'] = True
      if 'walltime' not in self.data:
        self.data['walltime']="00:72:00"
      if 'nodes' not in self.data:
        self.data['nodes']=1
      if 'ppn' not in self.data:
        self.data['ppn']=1
        

        
      nodes = self.data['nodes']
      ppn = self.data['ppn']
      walltime = self.data['walltime']
	

      self.data['addtlResc'][kw].append(value)      
      sf.write("#PBS -l walltime={}\n".format(walltime))
      sf.write("#PBS -l nodes={}:ppn={}\n".\
        format(str(nodes), str(ppn)))
      sf.write("#PBS -A {}\n".format(account))
      sf.write("#PBS -q {}\n".format(queue))
      
      print "queue type set automatically to PBS, set walltime and queue with %s.add_queue_resource(\"-q\", \"XXX\")" %(self.__name__)
    else:
      raise NotImplementedError
      
  def __setitem__(self, name, value):
    if name == 'queueType':
      if not (value in acceptableQueueTypes):
        raise ValueError("Undefined queueType %s, please\
        choose from this list:\n[%s]\n" %(value, ", ".join(acceptableQueueTypes)))
    self.data[name] = value
    
  def __getitem__(self, name):
    if not (name in self.data):
      print "Instance attribute '%s' is not defined" %(name)
      return
      
    return self.data[name]

  def __delitem__(self, name):
    del self.data[name]
  
  def add_wait(self, waitID):
    """ 
    Add passed job ID to list of jobs for this job submission to
    wait for. Can be called multiple times.
    """
    self.data['waitFor'].add(str(waitID))

  def add_Q_resource(self, kw, value):
    """
    Add passed keyword and value to a list of attributes that
    will be passed to the scheduler
    """
    self.data['addtlResc'][kw].append(value)

  def add_command(self, command):
    """
    add a single command-line executable
    input: string
    output: None
    """
    assert type(command) == str
    self.data['commandList'].append(command)

  def add_many_commands(self, commands):
    """
    add a many command-line executables
    input: iterable of strings
    output: None
    """    
    for command in commands:
      self.add_command(command)

  def show_commands(self):
    """
    print commands to stdout
    """
    for command in self.data['commandList']:
      print command
      
  def submit(self, command=None, **kwargs):
    """
    submit a job.
    if a command is given as the first argument, it will be appended to the list of args
    an shFile will be created if it doesn't exist
    """

    if command is not None:
      if type(command) == list:
        self.add_many_commands(command)
      elif type(command) == str:
        self.add_command(command)
      else:
        raise ValueError("I don't understand commands that aren't strings or lists")
      
    if len(self.data['commandList']) == 0:
      raise ValueError("I have no commands to run")
    
    if not (self.data['shWasWritten']):
      self.write_sh(**kwargs)
    
    if not os.path.exists(self.data['shFile']):
      raise UnboundLocalError("shFile is undefined")
  
    if self['queueType'] == "SGE" or self['queueType'] == "PBS" :
      p = subprocess.Popen(["qsub", self.data['shFile']], stdout=PIPE)
      output = p.communicate()[0].strip()
      try:
        jobId = re.findall(r'\d+', output)[0]
      except:
        raise RuntimeError("I didn't get back a well-formed job ID from the queue")
      self.data['sentJobs'].append(jobId)      
      return jobId
    else:
      raise UnboundLocalError("undefined queueType")

  def write_sh(self, queueType=None, useArray=None, chunks=1, 
               submit=False, shFile=None, jobName=None, joinArrayOut=True):
    """
    Create a file that is the shell script to be submitted to the
    job scheduler. 

    input kwargs:
    queueType: Scheduler type, either SGE or PBS    
    useArray: distribute jobs to array of jobs if True
    chunks: if array is True, Int number of commands per job
    submit: submit to scheduler default:False only write SH file if False
    shFile: name of shell file to write. must end with .sh
    jobName: name of this job
    joinArrayOut: set True if you want separate .e and .o files for each job

    returns:
    job ID assigned by scheduler for this job if submit is True,
    returns 0 if only writing SH file.
    
    """
   
      
    retVal = 0
    if shFile == None:
      if not ('shFile' in self.data):
        raise UnboundLocalError("shFile is undefined")
    else:
      self.data['shFile'] = shFile
      
    if not (self.data['shFile'].endswith(".sh")):
      raise ValueError("sh filename %s must end with .sh" %(shFile))


    if queueType == None:
      if not ('queueType' in self.data):
        raise UnboundLocalError("queueType is undefined")
    else:
      self.data['queueType'] = queueType


    if jobName == None:
      if not ('jobName' in self.data):
        raise UnboundLocalError("jobName is undefined")
    else:
      self.data['jobName'] = jobName

    if useArray is not None:
      self.data['useArray'] = useArray
    
    with open(self.data['shFile'], 'w') as sf: #memory-safe file opening

      if (self.data['queueType'] == 'SGE'):

        sf.write("#!/bin/sh\n")
        sf.write("##this script was generated by %s\n" %(os.path.abspath(__file__)))
        sf.write("#$ -N " + self.data['jobName']+"\n")

        sf.write("#$ -V\n")
        sf.write("#$ -S /bin/sh\n")
        sf.write("#$ -cwd\n")
      sf.write("#!/bin/sh\n")
      sf.write("#PBS -N "+self.data['jobName']+"\n")
      sf.write("#PBS -o "+self.data['jobName']+".out\n")
      sf.write("#PBS -e "+self.data['jobName']+".err\n")
      sf.write("#PBS -V\n")
      sf.write("cd $PBS_O_WORKDIR\n")

      if 'wait_for' in self.data:
        if self.data['wait_for']:   
          sf.write("#$ -hold_jid "+" ".join(self.data['wait_for'])+"\n")

      if 'addtl_resc' in self.data:
        if self.data['addtl_resc']:
          for key in self.data['addtl_resc']:
            for value in self.data['addtl_resc'][key]:
              sf.write("#$ "+key+" "+value+"\n")


              
        if len(list(self.data['waitFor'])) > 0:
          sf.write("#$ -hold_jid " + " ".join(self.data['waitFor'])+"\n")

        for key, values in self.data['addtlResc'].items():
          for value in values:
            sf.write("#$ "+key+" "+value+"\n")

        if ('useArray' in self.data) and (self.data['useArray'] == True):
          if joinArrayOut:
            sf.write("#$ -o " + self.data['jobName']+".out\n")
            sf.write("#$ -e " + self.data['jobName']+".err\n")          

          numberJobs = math.ceil(len(self.data['commandList'])/float(chunks))
          if self.maxJobs:
            sf.write("#$ -tc %d\n" %(self.maxJobs))
          sf.write("#$ -t 1-%d\n" %(numberJobs))            

          for i, cmd in enumerate(self.data['commandList']):
            sf.write("cmd[%d]=\"%s\"\n" %(i+1, cmd))

          arrayTaskId = "$SGE_TASK_ID"

          sf.write("let \"stop = (%s)*%d\"\n" %(arrayTaskId, chunks))
                         
          sf.write("let \"start = (%s-1)*(%d) + 1\"\n" %(arrayTaskId, chunks))
                
          sf.write("for i in `seq $start $stop`\n")
          sf.write("do\n")
          #sf.write("\tsleep $SGE_TASK_ID\n")
          sf.write("\teval ${cmd[$i]};\n")
          sf.write("done;\n")

        else:
          sf.write("#$ -o " + self.data['jobName']+".out\n")
          sf.write("#$ -e " + self.data['jobName']+".err\n")          
          for command in self.data['commandList']:
            sf.write(str(command) + "\n")
            
    self.shWasWritten=True #this instance has written an sh file

    if submit:
      return self.submit()
    return retVal
      




