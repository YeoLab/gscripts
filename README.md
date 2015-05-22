[![DOI](https://zenodo.org/badge/6604/YeoLab/gscripts.png)](http://dx.doi.org/10.5281/zenodo.12229)
## The code provided herein is without any warranty. 

It is mostly development code, bugs and inefficiencies are likely, please report any problems to the gscripts issues page.

## How to install gscripts

```
cd path/to/gscripts; easy_install . 
```
(NOT pip)

## How to import `qtools`

```
from gscripts.qtools import Submitter
```

For example, to submit a single line of code to a compute job, you can do:

```python
from gscripts.qtools import Submitter

exon_seqs = '/projects/ps-yeolab/obotvinnik/miso_helpers/hg19/se_exon2.fasta'
mirna_seqs = '/projects/ps-yeolab/genomes/mirbase/release_21/human_mature_17bp.fa'
rnahybrid_results = '/projects/ps-yeolab/obotvinnik/miso_helpers/hg19/se_exon2_RNAhybrid_mirbase_human_mature_17bp.txt'
command = 'RNAhybrid -c -s 3utr_human -q {} -t {} > {}'.format(mirna_seqs, exon_seqs, rnahybrid_results)
sub = Submitter([command], 'RNAhybrid', walltime='120:00:00', write_and_submit=True, nodes=1, ppn=1)
```

Which will write a file called `RNAhybrid.sh` which has these contents:

```bash
#!/bin/bash
#PBS -N RNAhybrid
#PBS -o RNAhybrid.sh.out
#PBS -e RNAhybrid.sh.err
#PBS -V
#PBS -l walltime=120:00:00
#PBS -l nodes=1:ppn=1
#PBS -A yeo-group
#PBS -q home

# Go to the directory from which the script was called
cd $PBS_O_WORKDIR
RNAhybrid -c -s 3utr_human -q /projects/ps-yeolab/genomes/mirbase/release_21/human_mature_17bp.fa -t /projects/ps-yeolab/obotvinnik/miso_helpers/hg19/se_exon2.fasta > /projects/ps-yeolab/obotvinnik/miso_helpers/hg19/se_exon2_RNAhybrid_mirbase_human_mature_17bp.txt
```

## How to import Python version of `which`

```
from gscripts import which
```
