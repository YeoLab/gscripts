from pybedtools import BedTool
from Bio import SeqIO
import os.path

annotations = "/home/jbrought/ce10_genes.bed"
path1 = "/oasis/tscc/scratch/jbrought/mirpipe/"
path2 = str(path1) + "v5_output/"

WT_libs = ["B1_NoIndex_L007_R1.R04_N9",
           "B1_NoIndex_L007_R1.R08_N16",
           "B1_NoIndex_L007_R1.R09_N10",
           "IC10.R10_N2_1",
           "IC10.R11_N2_2",
           "IC10.R12_N2_3"]

mut_libs = ["B1_NoIndex_L007_R1.R04_N9",
            "B1_NoIndex_L007_R1.R08_N16",
            "B1_NoIndex_L007_R1.R10_L10",
            "IC10.R13_L1_1",
            "IC10.R14_L1_2",
            "IC10.R15_L1_3"]

center = ".polyATrim.adapterTrim.groom."
tail = ".rg.sorted.bam"

#libs = [WT_libs, mut_libs]
libs = ["WT_libs.", "n2853_libs."]
lib_dic = {"WT_libs.": WT_libs, "n2853_libs.": mut_libs}


input_file = str(path1) + "sub.mature.fa"
input_mirs = SeqIO.parse(open(input_file), 'fasta')

names = []

for fasta in input_mirs:
    mirID = str(fasta.id)
    names.append(mirID)

for mir in names:
    for i in libs:
        curr_lib = [g + str(center) + str(mir) + str(tail) for g in lib_dic[i]]
        real_lib = []
        for item in curr_lib:
            if os.path.exists(item) == True:
                real_lib.append(item)
            else:
                print(str(item) + " Does not exist")
        x = BedTool()
        y = x.multi_bam_coverage(bams=real_lib, bed=annotations, s=True)
        outname = str(path2) + str(i) + str(mir) + ".genes.intersect.bed"
        y.saveas(str(outname))



    


