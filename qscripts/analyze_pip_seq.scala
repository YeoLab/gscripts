package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.picard.{ ReorderSam, SortSam, AddOrReplaceReadGroups, MarkDuplicates }
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.yeo._

class analyze_pip_seq extends QScript {
	@Input(doc = "input fastq file - or a list of fastq files", fullName="input", shortName="I")
	var input: File = _

	@Argument(doc = "species (hg19...)")
	var star_genome_location: String = _

	@Argument(doc = "location of chr_sizes")
	var chr_sizes: String = _

	@Argument(doc = "adapter to trim")
	var adapter: List[String] = Nil

	@Argument(doc = "STAR genome location")
	var star_genome_location: String = _
	
  case class sortSam(inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam {
    override def shortDescription = "sortSam"

    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP

  }

  case class star(input: File, output: File, genome_location: String) extends STAR {
       this.inFastq = input
       this.outSam = output
       this.genome = genome_location
       this.multimapNMax = 1
  }
  
  case class genomeCoverageBed(input: File, outBed : File, cur_strand : String) extends GenomeCoverageBed {
	this.inBam = input
	this.genomeSize = chr_sizes
	this.bedGraph = outBed
	this.strand = cur_strand
	this.split = true
  }

  case class star(input: File, output: File) extends STAR {
       this.inFastq = input
       this.outSam = output
       this.genome = star_genome_location
  }

  def script() {

	val fileList: Seq[File] = QScriptUtils.createSeqFromFile(input)
	var trackHubFiles: List[File] = List()

	for (fastq_file: File <- fileList) {

	 val noPolyAFastq = swapExt(fastq_file, ".fastq", ".polyATrim.fastq")
	 val noPolyAReport = swapExt(noPolyAFastq, ".fastq", ".report")

	 val noAdapterFastq = swapExt(noPolyAFastq, ".fastq", ".adapterTrim.fastq")
	 val adapterReport = swapExt(noAdapterFastq, ".fastq", ".report")

	 val filteredFastq = swapExt(noAdapterFastq, ".fastq", ".rmRep.fastq")
	 val filterd_results = swapExt(filteredFastq, ".fastq", ".counts")

	 val samFile = swapExt(filteredFastq, ".fastq", ".sam")
	 val sortedBamFile = swapExt(samFile, ".sam", ".sorted.bam")
	 val NRFFile = swapExt(sortedBamFile, ".bam", ".NRF")

	 val bedGraphFilePos = swapExt(sortedBamFile, ".bed", ".pos.bg")
	 val bigWigFilePos = swapExt(bedGraphFilePos, ".bg", ".bw")

	 val bedGraphFileNeg = swapExt(sortedBamFile, ".bed", ".neg.bg")
	 val bedGraphFileNegInverted = swapExt(bedGraphFileNeg, "neg.bg", "neg.t.bg")
	 val bigWigFileNeg = swapExt(bedGraphFileNegInverted, ".t.bg", ".bw")

	 add(new FastQC(inFastq = fastq_file))

	 add(new Cutadapt(inFastq = fastq_file, outFastq = noAdapterFastq, 
			 report = adapterReport, 
			 anywhere = adapter ++ List("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), 
			 overlap = 5, length = 18, quality_cutoff = 10))

	 add(new FilterRepetativeRegions(inFastq = noAdapterFastq, filterd_results, filteredFastq))

	 add(new FastQC(filteredFastq))
	 add(new star(filteredFastq, samFile))

	 add(sortSam(samFile, sortedBamFile, SortOrder.coordinate))

	 add(new CalculateNRF(inBam = sortedBamFile, genomeSize = chr_sizes, outNRF = NRFFile))

	 add(new genomeCoverageBed(input = sortedBamFile, outBed = bedGraphFilePos, cur_strand = "+"))
	 add(new BedGraphToBigWig(bedGraphFilePos, chr_sizes, bigWigFilePos))

	 add(new genomeCoverageBed(input = sortedBamFile, outBed = bedGraphFileNeg, cur_strand = "-"))
	 add(new NegBedGraph(inBedGraph = bedGraphFileNeg, outBedGraph = bedGraphFileNegInverted))
	 add(new BedGraphToBigWig(bedGraphFileNegInverted, chr_sizes, bigWigFileNeg))
	}
 }
}