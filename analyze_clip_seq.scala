package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.picard.{ReorderSam, SortSam, AddOrReplaceReadGroups, MarkDuplicates}
import org.broadinstitute.sting.queue.extensions.gatk._


class AnalizeCLIPSeq extends QScript {
  // Script argunment
  @Input(doc="input file")
  var input: List[File] = Nil

  @Argument(doc="species (hg19...)")
  var species: String =  _

  @Argument(doc="location of chr_sizes")
  var chr_sizes: String = _
 
  @Argument(doc="adapter to trim")
  var adapter: List[String] = Nil

  @Argument(doc="RBP binding modality is premRNA")
  var premRNA: Boolean = false

  @Argument(doc="location to place trackhub (must have the rest of the track hub made before starting script)")
  var location: String = _

 class FilterRepetativeRegions(@Input inFastq: File, @Output outCounts: File, @Output outNoRep: File) extends CommandLineFunction {

	def commandLine = "bowtie -S -q -p 4 -e 100 -l 20 --un %s all_ref %s | grep -v \"@\" | perl /nas3/yeolab/Software/pipeflower/count_aligned_from_sam.pl > %s".format(outNoRep, inFastq, outCounts)

 }

 class FastQC(@Input inFastq: File) extends CommandLineFunction {

	def commandLine = "fastqc %s".format(inFastq)

 }

 class CalculateNRF(@Input inBam: File, @Output outNRF: File, @Argument genomeSize: String) extends CommandLineFunction {

	def commandLine = "python /nas3/gpratt/gscripts/calculate_NRF.py " + required("--bam", inBam) + required("--genome", genomeSize) + " > " + outNRF
   
 }

 class Cutadapt(@Input inFastq: File, @Output outFastq: File, @Output report: File, @Argument anywhere: List[String] = Nil, front: List[String] = Nil, @Argument overlap: Option[Int] = None, error_rate: Option[Double] = None, length: Option[Int] = None, quality_cutoff: Option[Int] = None) extends CommandLineFunction {
	//see argunments on cutadapt command line for more documentation

	def commandLine = "cutadapt -f fastq --match-read-wildcards --times 2" + optional("-e", error_rate) + optional("-O", overlap) + optional("--quality-cutoff", quality_cutoff) + optional("-m", length) + repeat("-b", anywhere) + repeat("-f", front) + required("-o", outFastq) + required(inFastq) + " > " +  report
	this.isIntermediate = true
 }

 class MapWithSTAR(@Input inFastq: File, @Output samFile: File, @Argument genome: String) extends CommandLineFunction{

	def commandLine = "/nas3/yeolab/Software/STAR/STAR_2.3.0e/STAR --runMode alignReads --runThreadN 4 --genomeDir /nas3/yeolab/Software/STAR/genomes/2.2.0/%s --genomeLoad LoadAndRemove --readFilesIn %s --outSAMunmapped Within --outFilterMultimapNmax 1 --outFileNamePrefix %s --outStd SAM > %s".format(genome, inFastq, samFile, samFile)
	this.isIntermediate = true
 }

 class sortBam extends CommandLineFunction {
	@Input(doc="Bam file to sort") 
	var inBam: File = _

	@Output(doc="Sorted bam file", required=false)
	var outBam: File = _

	def commandLine = "samtools sort %s %s".format(inBam, outBam)	

 }

  class genomeCoverageBed(@Input inBam: File, @Argument genomeSize: String, @Output bedGraph: File, @Argument strand: String = null) extends CommandLineFunction{
	//When it comes time to refactor use the @Input(doc, requiered=False) pattern...
	def commandLine = "genomeCoverageBed " + optional("-strand", strand) + required("-split") + required("-bg") + required("-ibam", inBam) + required("-g", genomeSize) +  " > " + bedGraph
	this.isIntermediate = true
 }
 
  case class sortSam (inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam {
    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP

  }

  case class markDuplicates(inBam: File, outBam: File, metrics_file: File, remove_duplicates: Boolean) extends MarkDuplicates {
	this.input = List(inBam)
	this.output = outBam
	this.metrics = metrics_file
	this.REMOVE_DUPLICATES = remove_duplicates
	this.isIntermediate = false
 } 
 class bedGraphToBigWig(@Input inBedGraph: File, @Argument genomeSize: String, @Output bigWig: File) extends CommandLineFunction{

	def commandLine = "bedGraphToBigWig %s %s %s".format(inBedGraph, genomeSize, bigWig)

 }

 class Clipper(@Input inBam: File, @Argument species: String, @Argument premRNA: Boolean, @Output outBed: File) extends CommandLineFunction{

 	def commandLine = "clipper " +
			   required("-b", inBam) + 
		           required("-s", species) +
			   required("-o", outBed) + 
			   conditional(premRNA, "--premRNA") +  
			   required("--bonferroni") + 
			   required("--superlocal")
 }

 class Clip_Analysis(@Input inBam: File, @Input inBed: File, @Argument species: String, @Output metrics: File) extends CommandLineFunction {
	
	def commandLine = "clip_analysis " + 
			   required("--clusters", inBed) + 
			   required("-s", species) + 
			   required("--bam", inBam) + 
			   required("--regions_location", "/nas3/lovci/projects/ucscBED/%s".format(species)) +  
			   required("--AS_Structure", "/nas3/yeolab/Genome/ensembl/AS_STRUCTURE/%sdata4".format(species)) + 
			   required("--genome_location", "/nas3/yeolab/Genome/ucsc/%s/chromosomes/all.fa".format(species)) +
			   required("--phastcons_location", "/nas3/yeolab/Conservation/phastCons/hg19_46way/placentalMammals/reformat/hg19_phastcons.bw") +
			   required("--motifs",  "AAAAA") +
			   required("--nrand", 3) +
			   required("--runPhast") +
			   required("--runMotif") +
			   required("--runHomer") +
			   required("--metrics", metrics) 

 }

 class MakeTrackHub(@Input bwFiles: List[File], @Argument location: String) extends CommandLineFunction {
	//@Input(doc="Bam file to sort") 
	//var bwFiles: List[File] = Nil 
	
	//@Argument(doc="Location") 
	//var location: String = _
	
	def commandLine = "python /nas3/gpratt/gscripts/make_trackhubs.py" + repeat(bwFiles) + required("--location", location)

 	
 }

 class NegBedGraph(@Input inBedGraph: File, @Output outBedGraph: File) extends CommandLineFunction {
     def commandLine = "python /nas3/gpratt/gscripts/negBedGraph.py " + required("--bg", inBedGraph) + " > " + required(outBedGraph)     
     this.isIntermediate = true
 }
 
 class RunIDR(@Input inBam: File, @Output outResult: File, @Argument premRNA: Boolean, @Argument species: String, @Argument genome: String) extends CommandLineFunction {
     def commandLine = "python /nas3/gpratt/projects/idr/perform_idr.py " +
			required("--bam", inBam) +
			required("--out", outResult) +
			conditional(premRNA, "--premRNA") +
			required("--species", species) +
			required("--genome", genome) 
 }

  def script() {

//    val fileList: Seq[File] = QScriptUtils.createSeqFromFile(input)
     var trackHubFiles : List[File] = List()
   
    for (fastq_file: File <- input) {
    
    val noPolyAFastq  = swapExt(fastq_file, ".fastq", ".polyATrim.fastq")    
    val noPolyAReport = swapExt(noPolyAFastq, ".fastq", ".report")

    val noAdapterFastq = swapExt(noPolyAFastq, ".fastq", ".adapterTrim.fastq")
    val adapterReport = swapExt(noAdapterFastq, ".fastq", ".report")

    val filteredFastq = swapExt(noAdapterFastq, ".fastq", ".rmRep.fastq")    
    val filterd_results = swapExt(filteredFastq, ".fastq", ".counts")

    val samFile = swapExt(filteredFastq, ".fastq", ".sam")
    val sortedBamFile = swapExt(samFile, ".sam", ".sorted.bam")
    val NRFFile = swapExt(sortedBamFile, ".bam", ".NRF")

    val rmDupedBamFile = swapExt(sortedBamFile, ".bam", ".rmDup.bam")
    val rmDupedMetricsFile = swapExt(rmDupedBamFile, ".bam", ".metrics")  
    
    val bedGraphFilePos = swapExt(rmDupedBamFile, ".bam", ".pos.bg")  
    val bigWigFilePos   = swapExt(bedGraphFilePos, ".bg", ".bw")

    val bedGraphFileNeg = swapExt(rmDupedBamFile, ".bam", ".neg.bg")
    val bedGraphFileNegInverted = swapExt(bedGraphFileNeg, "neg.bg", "neg.t.bg")  
    val bigWigFileNeg   = swapExt(bedGraphFileNeg, ".t.bg", ".bw")
       
    val clipper_output = swapExt(rmDupedBamFile, ".bam", ".peaks.bed")
    val clipper_output_metrics = swapExt(clipper_output, ".bed", ".metrics")

    val IDRResult = swapExt(rmDupedBamFile, "", ".IDR")

    //add bw files to list for printing out later
    trackHubFiles = trackHubFiles ++ List(bedGraphFileNegInverted, bigWigFilePos)
    //trackHubFiles = trackHubFiles ++ List(clipper_output)
    add(new FastQC(inFastq=fastq_file)) 
    
    //filters out adapter reads
    add(new Cutadapt(inFastq=fastq_file, outFastq=noAdapterFastq, report=adapterReport, anywhere=adapter ++ List("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), overlap=5, length=18, quality_cutoff=10))
    add(new FilterRepetativeRegions(inFastq=noAdapterFastq, filterd_results, filteredFastq))
    add(new FastQC(filteredFastq))    
    add(new MapWithSTAR(filteredFastq, samFile, species))
    add(sortSam(samFile, sortedBamFile, SortOrder.coordinate))

    add(new CalculateNRF(inBam=sortedBamFile, genomeSize=chr_sizes, outNRF=NRFFile)) 
    add(markDuplicates(sortedBamFile, rmDupedBamFile, rmDupedMetricsFile, true)) 
	
    add(new genomeCoverageBed(inBam=rmDupedBamFile, genomeSize=chr_sizes, bedGraph=bedGraphFilePos, strand="+"))
    add(new bedGraphToBigWig(bedGraphFilePos, chr_sizes, bigWigFilePos))
    
    add(new genomeCoverageBed(inBam=rmDupedBamFile, genomeSize=chr_sizes, bedGraph=bedGraphFileNeg, strand="-"))
    add(new NegBedGraph(inBedGraph=bedGraphFileNeg, outBedGraph=bedGraphFileNegInverted))   
    add(new bedGraphToBigWig(bedGraphFileNegInverted, chr_sizes, bigWigFileNeg))

    add(new Clipper(inBam=rmDupedBamFile, species=species, outBed=clipper_output, premRNA=premRNA))
    add(new Clip_Analysis(rmDupedBamFile, clipper_output, species, clipper_output_metrics))
    
    add(new RunIDR(inBam=rmDupedBamFile, species=species, genome=chr_sizes, outResult=IDRResult, premRNA=premRNA))

  }
    add(new MakeTrackHub(trackHubFiles, location))   	  
 }
}





