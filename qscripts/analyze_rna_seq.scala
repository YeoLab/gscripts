package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.picard.{ ReorderSam, SortSam, AddOrReplaceReadGroups, MarkDuplicates }
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.yeo._

class AnalizeCLIPSeq extends QScript {
  // Script argunment
  @Input(doc = "input file or txt file of input files")
  var input: File = _

  @Argument(doc = "species (hg19...)")
  var species: String = _

  @Argument(doc = "location of chr_sizes")
  var chr_sizes: String = _

  @Argument(doc = "adapter to trim")
  var adapter: List[String] = Nil

  @Argument(doc = "flipped", required = false)
  var flipped: String = _

  @Argument(doc = "location to place trackhub (must have the rest of the track hub made before starting script)")
  var location: String = _

  case class sortSam(inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam {
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

  case class genomeCoverageBed(input: File, outBed : File, cur_strand: String) extends GenomeCoverageBed {
        this.inBam = input
        this.genomeSize = chr_sizes
        this.bedGraph = outBed
        this.strand = cur_strand
        this.split = true
  }

  case class singleRPKM(input: File, output: File, s: String) extends SingleRPKM {
	this.inCount = input
	this.outRPKM = output
	this.species = s
  }

  case class countTags(input: File, output: File, s: String) extends CountTags {
	this.inBam = input
	this.outCount = output
	this.species = s
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

      val bedGraphFilePos = swapExt(sortedBamFile, ".bam", ".pos.bg")
      val bigWigFilePos = swapExt(bedGraphFilePos, ".bg", ".bw")

      val bedGraphFileNeg = swapExt(sortedBamFile, ".bam", ".neg.bg")
      val bedGraphFileNegInverted = swapExt(bedGraphFileNeg, "neg.bg", "neg.t.bg")
      val bigWigFileNeg = swapExt(bedGraphFileNeg, ".t.bg", ".bw")

      val countFile = swapExt(sortedBamFile, "bam", "count")
      val RPKMFile = swapExt(countFile, "count", "RPKM")

      //add bw files to list for printing out later
      trackHubFiles = trackHubFiles ++ List(bedGraphFileNegInverted, bigWigFilePos)
      //trackHubFiles = trackHubFiles ++ List(clipper_output)
      add(new FastQC(inFastq = fastq_file))

      //filters out adapter reads
      add(new Cutadapt(inFastq = fastq_file, outFastq = noAdapterFastq, 
          report = adapterReport, 
          anywhere = adapter ++ List("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
        		  					 "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), 
          overlap = 5, length = 18, quality_cutoff = 6))
          
      add(new FilterRepetativeRegions(inFastq = noAdapterFastq, filterd_results, filteredFastq))
      add(new FastQC(filteredFastq))
      add(new STAR(filteredFastq, samFile, species))
      add(new sortSam(samFile, sortedBamFile, SortOrder.coordinate))

      add(new CalculateNRF(inBam = sortedBamFile, genomeSize = chr_sizes, outNRF = NRFFile))
      
      add(new genomeCoverageBed(input = sortedBamFile, outBed = bedGraphFilePos, cur_strand = "+"))
      add(new BedGraphToBigWig(bedGraphFilePos, chr_sizes, bigWigFilePos))

      add(new genomeCoverageBed(input = sortedBamFile, outBed = bedGraphFileNeg, cur_strand = "-"))
      add(new NegBedGraph(inBedGraph = bedGraphFileNeg, outBedGraph = bedGraphFileNegInverted))
      add(new BedGraphToBigWig(bedGraphFileNegInverted, chr_sizes, bigWigFileNeg))
	
      add(new countTags(input = sortedBamFile, output = countFile, s = species))	
      add(new singleRPKM(input = countFile, output = RPKMFile, s = species))

    }
    add(new MakeTrackHub(trackHubFiles, location))
  }
}





