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
  @Input(doc = "input file")
  var input: List[File] = Nil

  @Argument(doc = "species (hg19...)")
  var species: String = _

  @Argument(doc = "location of chr_sizes")
  var chr_sizes: String = _

  @Argument(doc = "adapter to trim")
  var adapter: List[String] = Nil

  @Argument(doc = "RBP binding modality is premRNA")
  var premRNA: Boolean = false

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

  def script() {

    //    val fileList: Seq[File] = QScriptUtils.createSeqFromFile(input)
    var trackHubFiles: List[File] = List()

    for (fastq_file: File <- input) {

      val noPolyAFastq = swapExt(fastq_file, ".fastq", ".polyATrim.fastq")
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
      val bigWigFilePos = swapExt(bedGraphFilePos, ".bg", ".bw")

      val bedGraphFileNeg = swapExt(rmDupedBamFile, ".bam", ".neg.bg")
      val bedGraphFileNegInverted = swapExt(bedGraphFileNeg, "neg.bg", "neg.t.bg")
      val bigWigFileNeg = swapExt(bedGraphFileNeg, ".t.bg", ".bw")

      val clipper_output = swapExt(rmDupedBamFile, ".bam", ".peaks.bed")
      val clipper_output_metrics = swapExt(clipper_output, ".bed", ".metrics")

      val IDRResult = swapExt(rmDupedBamFile, "", ".IDR")

      //add bw files to list for printing out later
      trackHubFiles = trackHubFiles ++ List(bedGraphFileNegInverted, bigWigFilePos)
      //trackHubFiles = trackHubFiles ++ List(clipper_output)
      add(new FastQC(inFastq = fastq_file))

      //filters out adapter reads
      add(new Cutadapt(inFastq = fastq_file, outFastq = noAdapterFastq, 
          report = adapterReport, 
          anywhere = adapter ++ List("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), 
          overlap = 5, length = 18, quality_cutoff = 10))
          
      add(new FilterRepetativeRegions(inFastq = noAdapterFastq, filterd_results, filteredFastq))
      add(new FastQC(filteredFastq))
      add(new STAR(filteredFastq, samFile, species))
      add(sortSam(samFile, sortedBamFile, SortOrder.coordinate))

      add(new CalculateNRF(inBam = sortedBamFile, genomeSize = chr_sizes, outNRF = NRFFile))
      add(markDuplicates(sortedBamFile, rmDupedBamFile, rmDupedMetricsFile, true))

      add(new GenomeCoverageBed(inBam = rmDupedBamFile, genomeSize = chr_sizes, bedGraph = bedGraphFileNeg, strand = "-"))
      add(new NegBedGraph(inBedGraph = bedGraphFileNeg, outBedGraph = bedGraphFileNegInverted))
      add(new BedGraphToBigWig(bedGraphFileNegInverted, chr_sizes, bigWigFileNeg))

      add(new Clipper(inBam = rmDupedBamFile, species = species, outBed = clipper_output, premRNA = premRNA))
      add(new ClipAnalysis(rmDupedBamFile, clipper_output, species, clipper_output_metrics))

      add(new IDR(inBam = rmDupedBamFile, species = species, genome = chr_sizes, outResult = IDRResult, premRNA = premRNA))

    }
    add(new MakeTrackHub(trackHubFiles, location))
  }
}





