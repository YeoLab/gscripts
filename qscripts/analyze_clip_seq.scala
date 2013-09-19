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

//  @Argument(doc = "species (hg19...)")
//  var species: String = _
  
//  @Argument(doc = "STAR genome location")
//  var star_genome_location: String = _

//  @Argument(doc = "location of chr_sizes")
//  var chr_sizes: String = _

  @Argument(doc = "adapter to trim")
  var adapter: List[String] = Nil

  @Argument(doc = "RBP binding modality is premRNA")
  var premRNA: Boolean = false

  @Argument(doc = "read ids have randomers")
  var barcoded: Boolean = false

  @Argument(doc = "location to place trackhub (must have the rest of the track hub made before starting script)", required=false)
  var location: String = "clip_seq"

  case class sortSam(inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam {
    override def shortDescription = "sortSam"
    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP

  }


  class RemoveDuplicates(@Input inBam: File, @Output outResult: File, @Argument metrics_file: String) extends CommandLineFunction {
  override def shortDescription = "RemoveDuplicates"
  def commandLine = "python ~/gscripts/gscripts/clipseq/barcode_collapse.py " +
    required("--bam", inBam) +
    required("--out_file", outResult) +
    conditional(barcoded, "--randomer") +
    required("--metrics_file", metrics_file) 
  }

  class FixScores(@Input inBed: File, @Output outBed: File) extends CommandLineFunction {
  override def shortDescription = "FixScores"
  def commandLine = "python ~/gscripts/gscripts/clipseq/fix_scores.py " +
    required("--bed", inBed) +
    required("--out_file", outBed)
    this.isIntermediate = true
  }

  class Ripseeker(@Input inBam: File, @Output outBed: File) extends CommandLineFunction {
   override def shortDescription = "ripseeker"
   def commandLine = "python ~/gscripts/gscripts/clipseq/run_ripseeker.py " + 
       		     required("--bam", inBam) +
   		     required("--out", outBed) 
  }
   
  class BamToBed(@Input inBam: File, @Output outBed: File) extends CommandLineFunction {
  override def shortDescription = "bamToBed"
   def commandLine = "bamToBed " +
       		     required("-i", inBam) + 
		     required("-split") + " > " + outBed
   this.isIntermediate = true

  }

  class Pyicoclip(@Input inBed: File, @Output outBed: File, @Argument regions: String) extends CommandLineFunction {
       override def shortDescription = "pyicolip"
       def commandLine =  "pyicoclip " +
       	   	       	  required(inBed) + 
			  required(outBed) +
			  required("-f",  "bed") + 
			  required("--region", regions)
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

  case class genomeCoverageBed(input: File, outBed : File, cur_strand : String, genome : String) extends GenomeCoverageBed {
        this.inBam = input
        this.genomeSize = genome
        this.bedGraph = outBed
        this.strand = cur_strand
        this.split = true
  }

  case class star(input: File, output: File, genome_location: String) extends STAR {
       this.inFastq = input
       this.outSam = output
       this.genome = genome_location
       this.multimapNMax = 1
  }

  def starGenomeLocation(genome: String) : String = {
  //Returns star genome Location for TSCC, could eventually be factored out into conf file
   
   var retval = "none"
   if (genome == "hg19") {
      retval = "/projects/ps-yeolab/genomes/hg19/star" 
   }else if(genome == "mm9") {
      retval = "/projects/ps-yeolab/genomes/mm9/star"
   }
   
   retval
  }

  def chromSizeLocation(genome: String) : String = {
  //Returns star genome Location for TSCC, could eventually be factored out into conf file
   
   var retval = "none"
   if (genome == "hg19") {
      retval = "/projects/ps-yeolab/genomes/hg19/hg19.chrom.sizes" 
   }else if(genome == "mm9") {
      retval = "/projects/ps-yeolab/genomes/mm9/mm9.chrom.sizes"
   }
   
   retval
  }

  def regionsLocation(genome: String) : String = {
  //Returns star genome Location for TSCC, could eventually be factored out into conf file
   
   var retval = "none"
   if (genome == "hg19") {
      retval = "/projects/ps-yeolab/genomes/hg19/hg19data4" 
   }else if(genome == "mm9") {
      retval = "/projects/ps-yeolab/genomes/mm9/mm9data4"
   }
   
   retval
  }

  def asStructureLocation(genome: String) : String = {
  //Returns star genome Location for TSCC, could eventually be factored out into conf file
   
   var retval = "none"
   if (genome == "hg19") {
      retval = "/projects/ps-yeolab/genomes/hg19/hg19data4" 
   }else if(genome == "mm9") {
      retval = "/projects/ps-yeolab/genomes/mm9/mm9data4"
   }
   
   retval
  }

  def genomeLocation(genome: String) : String = {
  //Returns star genome Location for TSCC, could eventually be factored out into conf file
   
   var retval = "none"
   if (genome == "hg19") {
      retval = "/projects/ps-yeolab/genomes/hg19/chromosomes/all.fa" 
   }else if(genome == "mm9") {
      retval = "/projects/ps-yeolab/genomes/mm9/chromosomes/all.fa" 
   }
   
   retval
  }

  def phastconsLocation(genome: String) : String = {
  //Returns star genome Location for TSCC, could eventually be factored out into conf file
   
   var retval = "none"
   if (genome == "hg19") {
      retval = "/projects/ps-yeolab/genomes/hg19/hg19_phastcons.bw" 
   }else if(genome == "mm9") {
      retval = "/projects/ps-yeolab/genomes/mm9/mm9_phastcons.bw" 
   }
   
   retval
  }

  def genicRegionsLocation(genome: String) : String = {
  //Returns star genome Location for TSCC, could eventually be factored out into conf file
   
   var retval = "none"
   if (genome == "hg19") {
      retval = "/home/gpratt/clipper/clipper/data/hg19.AS.STRUCTURE.COMPILED.bed" 
   }else if(genome == "mm9") {
      retval = "/home/gpratt/clipper/clipper/data/mm9.AS.STRUCTURE.COMPILED.bed" 
   }
   
   retval
  }

  def gffDbLocation(genome: String) : String = {
  //Returns star genome Location for TSCC, could eventually be factored out into conf file
   
   var retval = "none"
   if (genome == "hg19") {
      retval = "/projects/ps-yeolab/genomes/hg19/gencode.v17.annotation.gtf.db" 
   }else if(genome == "mm9") {
      retval = "/projects/ps-yeolab/genomes/mm9/Mus_musculus.NCBIM37.64.fixed.gtf.db" 
   }
   
   retval
  }
  

  def script() {

    val fileList = QScriptUtils.createArgsFromFile(input)
    var trackHubFiles: List[File] = List()

    for (item : Tuple2[File, String] <- fileList) {
      var fastq_file: File = item._1
      var genome: String = item._2 
 
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
      val bigWigFileNeg = swapExt(bedGraphFileNeg, ".bg", ".norm.bw")
      val bedGraphFileNegInverted = swapExt(bedGraphFileNeg, "neg.bg", "neg.t.bg")
      val bigWigFileNegInverted = swapExt(bedGraphFileNegInverted, ".t.bg", ".bw")

      val clipper_output = swapExt(rmDupedBamFile, ".bam", ".peaks.bed")
      val fixed_clipper_output = swapExt(clipper_output, ".bed", ".fixed.bed")
      val bigBed_output = swapExt(fixed_clipper_output, ".bed", ".bb")
      val clipper_output_metrics = swapExt(clipper_output, ".bed", ".metrics")
      
      val rmDupedBedFile = swapExt(rmDupedBamFile, ".bam", ".bed")
      val pyicoclipResults = swapExt(rmDupedBedFile, ".bed", ".pyicoclip.bed")
      val ripseekerResults = swapExt(rmDupedBamFile, ".bam", ".ripseeker.bed")
      val IDRResult = swapExt(rmDupedBamFile, "", ".IDR")

      //add bw files to list for printing out later
      trackHubFiles = trackHubFiles ++ List(bedGraphFileNegInverted, bigWigFilePos)
      trackHubFiles = trackHubFiles ++ List(bigBed_output)
//      add(new FastQC(inFastq = fastq_file))

      //filters out adapter reads
      add(new Cutadapt(inFastq = fastq_file, outFastq = noAdapterFastq, 
          report = adapterReport, 
          anywhere = adapter ++ List("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
        		  	     "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"), 
          overlap = 5, length = 18, quality_cutoff = 6))
      
      

      add(new FilterRepetativeRegions(inFastq = noAdapterFastq, filterd_results, filteredFastq))
      add(new FastQC(filteredFastq))
      add(star(input = filteredFastq, output = samFile, genome_location = starGenomeLocation(genome)))
      add(sortSam(samFile, sortedBamFile, SortOrder.coordinate))
      add(addOrReplaceReadGroups(sortedBamFile, rgSortedBamFile))

      add(new CalculateNRF(inBam = rgSortedBamFile, genomeSize = chromSizeLocation(genome), outNRF = NRFFile))
      add(new RemoveDuplicates(rgSortedBamFile, rmDupedBamFile, rmDupedMetricsFile))

      add(new genomeCoverageBed(input = rmDupedBamFile, outBed = bedGraphFilePos, cur_strand = "+", genome = chromSizeLocation(genome)))
      add(new BedGraphToBigWig(bedGraphFilePos, chromSizeLocation(genome), bigWigFilePos))

      add(new genomeCoverageBed(input = rmDupedBamFile, outBed = bedGraphFileNeg, cur_strand = "-", genome = chromSizeLocation(genome)))
      add(new BedGraphToBigWig(bedGraphFileNeg, chromSizeLocation(genome), bigWigFileNeg))
      add(new NegBedGraph(inBedGraph = bedGraphFileNeg, outBedGraph = bedGraphFileNegInverted))
      add(new BedGraphToBigWig(bedGraphFileNegInverted, chromSizeLocation(genome), bigWigFileNegInverted))

      add(new Clipper(inBam = rmDupedBamFile, species = genome, outBed = clipper_output, premRNA = premRNA))
            
      add(new FixScores(inBed = clipper_output, outBed = fixed_clipper_output))
      
      add(new BedToBigBed(inBed = fixed_clipper_output, genomeSize = chromSizeLocation(genome), outBigBed = bigBed_output))

      add(new ClipAnalysis(rmDupedBamFile, clipper_output, genome, clipper_output_metrics, regions_location = regionsLocation(genome),
     	      		   AS_Structure = asStructureLocation(genome), genome_location = genomeLocation(genome), 
			   phastcons_location = phastconsLocation(genome), gff_db = gffDbLocation(genome),
			   bw_pos=bigWigFilePos, bw_neg=bigWigFileNeg))

      add(new BamToBed(inBam=rmDupedBamFile, outBed=rmDupedBedFile))
      add(new Pyicoclip(inBed = rmDupedBedFile, outBed = pyicoclipResults, regions = genicRegionsLocation(genome) ))
      add(new Ripseeker(inBam=rmDupedBamFile, outBed=ripseekerResults))
      add(new IDR(inBam = rmDupedBamFile, species = genome, genome = chromSizeLocation(genome), outResult = IDRResult, premRNA = premRNA))

    }
//    add(new MakeTrackHub(trackHubFiles, location, genome))
  }
}