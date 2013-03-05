package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.picard.{ReorderSam, SortSam, AddOrReplaceReadGroups}
import org.broadinstitute.sting.queue.extensions.gatk._


class HelloWorld extends QScript {
  // Script argunment
  @Input(doc="input file")
  var input: File = _

  @Argument(doc="species (hg19...)")
  var species: String =  _

  @Argument(doc="location of chr_sizes")
  var chr_sizes: String = _
 
 class FilterRepetativeRegions(@Input inFasta: File, @Output outCounts: File, @Output outNoRep: File) extends CommandLineFunction {

	def commandLine = "bowtie -S -q -p 4 -e 100 -l 20 --un %s all_ref %s | grep -v \"@\" | perl /nas3/yeolab/Software/pipeflower/count_aligned_from_sam.pl > %s".format(outNoRep, inFasta, outCounts)

 }

 class FastQC(@Input inFasta: File) extends CommandLineFunction {

	def commandLine = "fastqc %s".format(inFasta)

 }

 class cutadapt(@Input inFasta: File, @Output outFasta: File, @Output report: File, @Argument anywhere: List[String], @Argument front: List[String], @Argument overlap: Int, @Argument error_rate: Float) extends CommandLineFunction {
	#see argunments on cutadapt command line for more documentation
	def commandLine = "cutadapt -f fastq" + optional("-e", error_rate) + optinal("-O", overlap) + repeat("-b", anywhere) + repeat("-f", front) + required("-o", outFasta) + required(">", report)
 }

 class MapWithSTAR(@Input inFasta: File, @Output samFile: File, @Argument genome: String) extends CommandLineFunction{

	def commandLine = "/nas3/yeolab/Software/STAR/STAR_2.3.0e/STAR --runMode alignReads --runThreadN 4 --genomeDir /nas3/yeolab/Software/STAR/genomes/2.2.0/%s --genomeLoad LoadAndRemove --readFilesIn %s --outSAMunmapped Within --outFilterMultimapNmax 10 --outStd SAM > %s".format(genome, inFasta, samFile)

 }

 class sortBam extends CommandLineFunction {
	@Input(doc="Bam file to sort") 
	var inBam: File = _

	@Output(doc="Sorted bam file", required=false)
	var outBam: File = _

	def commandLine = "samtools sort %s %s".format(inBam, outBam)	

 }
  class genomeCoverageBed(@Input inBam: File, @Argument genomeSize: String, @Output bedGraph: File) extends CommandLineFunction{

	def commandLine = "genomeCoverageBed -split -bg -ibam %s -g %s > %s".format(inBam, genomeSize, bedGraph)

}
 
  case class sortSam (inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam {
    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP

  }


 class bedGraphToBigWig(@Input inBedGraph: File, @Argument genomeSize: String, @Output bigWig: File) extends CommandLineFunction{

	def commandLine = "bedGraphToBigWig %s %s %s".format(inBedGraph, genomeSize, bigWig)

}

 class Clipper(@Input inBam: File, @Argument species: String, @Output outBed: File) extends CommandLineFunction{

 	def commandLine = "clipper -b %s -s %s -o %s".format(inBam, species, outBed) 

 }

 class Clip_Analysis(@Input inBam: File, @Input inBed: File, @Argument species: String) extends CommandLineFunction {
	
	def commandLine = "python /nas3/gpratt/clipper/clipper/src/CLIP_analysis.py --clusters %s -s %s --bam %s --regions_location /nas3/lovci/projects/ucscBED/%s --AS_Structure /nas3/yeolab/Genome/ensembl/AS_STRUCTURE/%sdata4 --genome_location /nas3/yeolab/Genome/ucsc/%s/chromosomes/all.fa --motif AAAAA --nrand 3 --rePhast".format(inBed, species, inBam, species, species, species) 

 }

  def script() {

    val fileList: Seq[File] = QScriptUtils.createSeqFromFile(input)

    for (fasta_file: File <- fileList) {
    
    val noPolyAFastq  = swapExt(fastq_file, ".fastq", ".polyATrim.fastq")
    val filteredFastq = swapExt(noPolyAFastq, ".fastq", ".rmRep.fastq")
    val noPolyAReport = swapExt(noPolyAFastq, ".fastq", ".report")
    val filterd_results = swapExt(fastq_file, ".fastq", ".counts")
    val samFile = swapExt(filteredFasta, ".fastq", ".sam")
    val sortedBamFile = swapExt(samFile, ".sam", ".sorted.bam")  
    val bedGraphFile = swapExt(sortedBamFile, ".bam", ".bg")  
    val bigWigFile   = swapExt(bedGraphFile, ".bg", ".bw")
    val clipper_output = swapExt(sortedBamFile, ".bam", ".peaks.bed")
 
    add(new FastQC(inFasta=fasta_file, outFasta=noPolyAFastq, report=noPolyAReport, anywhere=List("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT", overlap=7, error_rate=0.03))
    #filters out poly-a tails (and maybe other types of poly in the future)
    add(new Cutadapt(fasta_file, 
    add(new FilterRepetativeRegions(fasta_file, filterd_results, filteredFasta))
    add(new FastQC(filteredFasta))    
    add(new MapWithSTAR(filteredFasta, samFile, species))
    add(sortSam(samFile, sortedBamFile, SortOrder.coordinate))
    add(new genomeCoverageBed(sortedBamFile, chr_sizes, bedGraphFile))
    add(new bedGraphToBigWig(bedGraphFile, chr_sizes, bigWigFile))

    add(new Clipper(sortedBamFile, species, clipper_output))
    add(new Clip_Analysis(sortedBamFile, clipper_output, species))
} 	  
 }
}





