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
  var input: File = _

  @Argument(doc = "species (hg19...)")
  var species: String = _
  
  @Argument(doc = "STAR genome location")
  var star_genome_location: String = _

  @Argument(doc = "location of chr_sizes")
  var chr_sizes: String = _

  @Argument(doc = "adapter to trim")
  var adapter: List[String] = Nil

  @Argument(doc = "RBP binding modality is premRNA")
  var premRNA: Boolean = false

  @Argument(doc = "location to place trackhub (must have the rest of the track hub made before starting script)", required=false)
  var location: String = "clip_seq"

  @Argument(doc = "location of genomic regions (AS_Structure proxintron... generally in /nas3/yeolab/Genome/ensembl/AS_STRUCTURE/xxdata4)")
  var regions_location: String = _

  @Argument(doc = "AS STRUCTURE Files")
  var AS_Structure: String = _

  @Argument(doc = ".fa file for full genome")
  var genome_location: String = _
 
  @Argument(doc = "location of phastcons conservation bw file")
  var phastcons_location: String = _
  
  case class sortSam(inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam {
    override def shortDescription = "sortSam"
    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP

  }

  case class markDuplicates(inBam: File, outBam: File, metrics_file: File) extends MarkDuplicates {
    override def shortDescription = "MarkDuplicates"
    this.input = List(inBam)
    this.output = outBam
    this.metrics = metrics_file
    this.REMOVE_DUPLICATES = true
    this.isIntermediate = false
  }
  
  case class addOrReplaceReadGroups(inBam : File, outBam : File) extends AddOrReplaceReadGroups {
       override def shortDescription = "AddOrReplaceReadGroups"
       this.input = List(inBam)
       this.output = outBam
       this.RGLB = "foo" //should record library id
       this.RGPL = "illumina"
       this.RGPU = "bar" //can record barcode information (once I get bardoding figured out)
       this.RGSM = "baz"
       this.RGCN = "Biogem" //need to add switch between biogem and singapore
       this.RGID = "1"
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
       this.multimapNMax = 1
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
      val rgSortedBamFile = swapExt(sortedBamFile, ".sorted.bam", ".sorted.rg.bam")

      val NRFFile = swapExt(rgSortedBamFile, ".bam", ".NRF")

      val rmDupedBamFile = swapExt(rgSortedBamFile, ".bam", ".rmDup.bam")
      val rmDupedMetricsFile = swapExt(rmDupedBamFile, ".bam", ".metrics")

      val bedGraphFilePos = swapExt(rmDupedBamFile, ".bam", ".pos.bg")
      val bigWigFilePos = swapExt(bedGraphFilePos, ".bg", ".bw")

      val bedGraphFileNeg = swapExt(rmDupedBamFile, ".bam", ".neg.bg")
      val bedGraphFileNegInverted = swapExt(bedGraphFileNeg, "neg.bg", "neg.t.bg")
      val bigWigFileNeg = swapExt(bedGraphFileNegInverted, ".t.bg", ".bw")

      val clipper_output = swapExt(rmDupedBamFile, ".bam", ".peaks.bed")
      val bigBed_output = swapExt(clipper_output, ".bed", ".bb")
      val clipper_output_metrics = swapExt(clipper_output, ".bed", ".metrics")

      val IDRResult = swapExt(rmDupedBamFile, "", ".IDR")

      //add bw files to list for printing out later
      trackHubFiles = trackHubFiles ++ List(bedGraphFileNegInverted, bigWigFilePos)
      trackHubFiles = trackHubFiles ++ List(bigBed_output)
      add(new FastQC(inFastq = fastq_file))

      //filters out adapter reads
      add(new Cutadapt(inFastq = fastq_file, outFastq = noAdapterFastq, 
          report = adapterReport, 
          anywhere = adapter ++ List("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
        		  					 "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), 
          overlap = 5, length = 18, quality_cutoff = 6))
      
      

      add(new FilterRepetativeRegions(inFastq = noAdapterFastq, filterd_results, filteredFastq))
      add(new FastQC(filteredFastq))
      add(star(input = filteredFastq, output = samFile))
      add(sortSam(samFile, sortedBamFile, SortOrder.coordinate))
      add(addOrReplaceReadGroups(sortedBamFile, rgSortedBamFile))

      add(new CalculateNRF(inBam = rgSortedBamFile, genomeSize = chr_sizes, outNRF = NRFFile))
      add(markDuplicates(rgSortedBamFile, rmDupedBamFile, rmDupedMetricsFile))

      add(new genomeCoverageBed(input = rmDupedBamFile, outBed = bedGraphFilePos, cur_strand = "+"))
      add(new BedGraphToBigWig(bedGraphFilePos, chr_sizes, bigWigFilePos))

      add(new genomeCoverageBed(input = rmDupedBamFile, outBed = bedGraphFileNeg, cur_strand = "-"))
      add(new NegBedGraph(inBedGraph = bedGraphFileNeg, outBedGraph = bedGraphFileNegInverted))
      add(new BedGraphToBigWig(bedGraphFileNegInverted, chr_sizes, bigWigFileNeg))

      add(new Clipper(inBam = rmDupedBamFile, species = species, outBed = clipper_output, premRNA = premRNA))

      add(new BedToBigBed(inBed = clipper_output, genomeSize = chr_sizes, outBigBed = bigBed_output))

      add(new ClipAnalysis(rmDupedBamFile, clipper_output, species, clipper_output_metrics, regions_location = regions_location,
      	      		   AS_Structure = AS_Structure, genome_location = genome_location, phastcons_location = phastcons_location))

      add(new IDR(inBam = rmDupedBamFile, species = species, genome = chr_sizes, outResult = IDRResult, premRNA = premRNA))

    }
    add(new MakeTrackHub(trackHubFiles, location, species))
  }
}





