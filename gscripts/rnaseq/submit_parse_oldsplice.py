import glob

files = glob.glob("*.splices")

samples = []
filenames = []

species = None
try:
    species = sys.argv[1]

except:
    print "usage: submit_parse_oldsplice.py <species>"


assert len(species) > 0

for filename in files:

    filenames.append(filename)
    sample = filename.replace(".splices", "").replace(".flip", "_flip")
    samples.append(sample)


from gscripts.qtools import Submitter

sub = Submitter()

cmd = "parse_oldsplice.py --species %s" %species

for filename, sample in zip(filenames, samples):
    cmd += " --sample %s %s " %(filename, sample)

print cmd


cmd = [cmd]


sub.job(command_list=cmd, array=False, sh_file="parse.sh", job_name="parse", submit=False, queue="home", ppn=1)
