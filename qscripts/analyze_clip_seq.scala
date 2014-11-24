package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.picard.{ ReorderSam, SortSam, AddOrReplaceReadGroups, MarkDuplicates }
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.yeo._
import org.broadinstitute.sting.queue.util.TsccUtils._
class AnalizeCLIPSeq extends QScript {
  // Script argunment
  @Input(doc = "input file")
  var input: File = _

  @Argument(doc = "adapter to trim")
  var adapter: List[String] = Nil

  @Argument(doc = "RBP binding modality is premRNA")
  var premRNA: Boolean = false

  @Argument(doc = "reads are on opposite strand")
  var reverse_strand: Boolean = false

  @Argument(doc = "read ids have randomers")
  var barcoded: Boolean = false

  @Argument(doc = "use alpha version of STAR")
  var alpha: Boolean = false

  @Argument(doc = "start processing from uncollapsed bam file")
  var fromBam: Boolean = false

  case class clipper(in: File, out: File, genome: String, isPremRNA: Boolean, reverse: Boolean) extends Clipper 
  {
    
    
    this.inBam = in
    this.outBed = out
    this.species = genome
    this.premRNA = isPremRNA
    this.superlocal = superlocal
    this.reverse_strand = reverse
  }

  case class cutadapt(fastq_file: File, noAdapterFastq: File, adapterReport: File, adapter: List[String]) extends Cutadapt{
       override def shortDescription = "cutadapt"

       this.inFastq = fastq_file
       this.outFastq = noAdapterFastq 
       this.report = adapterReport
       this.anywhere = adapter ++ List("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 
        		  	     "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT") 
       this.overlap = 5
       this.length = 18
       this.quality_cutoff = 6 
       this.isIntermediate = true
  }

  case class filterRepetitiveRegions(noAdapterFastq: File, filteredResults: File, filteredFastq: File) extends FilterRepetitiveRegions {
       override def shortDescription = "FilterRepetitiveRegions"
       
       this.inFastq = noAdapterFastq
       this.outCounts = filteredResults
       this.outNoRep = filteredFastq
       this.isIntermediate = true
  }

  case class sortSam(inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam {
    override def shortDescription = "sortSam"
    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP
    this.isIntermediate = false
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
       this.isIntermediate = true
  } 

  case class genomeCoverageBed(input: File, outBed : File, cur_strand : String, genome : String) extends GenomeCoverageBed {
        this.inBam = input
        this.genomeSize = genome
        this.bedGraph = outBed
        this.strand = cur_strand
        this.split = true
	this.isIntermediate = true
  }

  case class star(input: File, output: File, genome_location: String) extends STAR {
       this.inFastq = input
       this.outSam = output
       this.genome = genome_location
       this.multimapNMax = 1
       this.isIntermediate = true
       this.alpha = alpha
  }

  case class samtoolsIndexFunction(input: File, output: File) extends SamtoolsIndexFunction {
       override def shortDescription = "indexBam"
       this.bamFile = input
       this.bamFileIndex = output
  }

  case class samtoolsMergeFunction(inBams: Seq[File], outBam: File) extends SamtoolsMergeFunction {
       override def shortDescription = "samtoolsMerge"
       this.inputBams = inBams
       this.outputBam = outBam
  }

  case class singleRPKM(input: File, output: File, s: String) extends SingleRPKM {
        this.inCount = input
        this.outRPKM = output
  }

  case class countTags(input: File, index: File, output: File, a: String) extends CountTags {
        this.inBam = input
        this.outCount = output
        this.tags_annotation = a
  }

// performs downstream analysis on a bam file of interest, peak calling, generating RPKMs and all other types of analysis  
  def downstream_analysis(bamFile : File, bamIndex: File, genome : String) {

      	     val bedGraphFilePos = swapExt(bamFile, ".bam", ".pos.bg")
	     val bedGraphFilePosNorm = swapExt(bedGraphFilePos, ".pos.bg", ".norm.pos.bg")
      	     val bigWigFilePos = swapExt(bedGraphFilePosNorm, ".bg", ".bw")

      	     val bedGraphFileNeg = swapExt(bamFile, ".bam", ".neg.bg")
	     val bedGraphFileNegNorm = swapExt(bedGraphFileNeg, ".neg.bg", ".norm.neg.bg")
      	     val bedGraphFileNegInverted = swapExt(bedGraphFileNegNorm, "neg.bg", "neg.t.bg")
      	     val bigWigFileNegInverted = swapExt(bedGraphFileNegInverted, ".t.bg", ".bw")

      	     val clipper_output = swapExt(bamFile, ".bam", ".peaks.bed")
	     val kasey_output = swapExt(bamFile, ".bam", ".peaks.kasey.bed")
      	     val fixed_clipper_output = swapExt(clipper_output, ".bed", ".fixed.bed")
      	     val bigBed_output = swapExt(fixed_clipper_output, ".bed", ".bb")
      	     val clipper_output_metrics = swapExt(clipper_output, ".bed", ".metrics")
	     val kasey_output_metrics = swapExt(kasey_output, ".bed", ".metrics")
	     
	     val rmDupedBedFile = swapExt(bamFile, ".bam", ".bed")
      	     val pyicoclipResults = swapExt(bamFile, ".bed", ".pyicoclip.bed")
      	     val IDRResult = swapExt(bamFile, "", ".IDR")
	     
      	     val countFile = swapExt(bamFile, "bam", "count")
      	     val RPKMFile = swapExt(countFile, "count", "RPKM")

	     //add bw files to list for printing out later

      	     add(new genomeCoverageBed(input = bamFile, outBed = bedGraphFilePos, cur_strand = "+", genome = chromSizeLocation(genome)))
	     add(new NormalizeBedGraph(inBedGraph = bedGraphFilePos, inBam = bamFile, outBedGraph = bedGraphFilePosNorm))
      	     add(new BedGraphToBigWig(bedGraphFilePosNorm, chromSizeLocation(genome), bigWigFilePos))

      	     add(new genomeCoverageBed(input = bamFile, outBed = bedGraphFileNeg, cur_strand = "-", genome = chromSizeLocation(genome)))
      	     add(new NormalizeBedGraph(inBedGraph = bedGraphFileNeg, inBam = bamFile, outBedGraph = bedGraphFileNegNorm))
	     add(new NegBedGraph(inBedGraph = bedGraphFileNegNorm, outBedGraph = bedGraphFileNegInverted))
      	     add(new BedGraphToBigWig(bedGraphFileNegInverted, chromSizeLocation(genome), bigWigFileNegInverted))
	     
      	     add(new clipper(in = bamFile, genome = genome, out = clipper_output, isPremRNA = premRNA, reverse=reverse_strand))
	     
      	     add(new FixScores(inBed = clipper_output, outBed = fixed_clipper_output))

	     add(new BedToBigBed(inBed = fixed_clipper_output, genomeSize = chromSizeLocation(genome), outBigBed = bigBed_output))

      	     add(new ClipAnalysis(bamFile, clipper_output, genome, clipper_output_metrics, AS_Structure = asStructureLocation(genome), 
	     			  genome_location = genomeLocation(genome), phastcons_location = phastconsLocation(genome), 
	     			  gff_db = gffDbLocation(genome)))

	     add(new BamToBed(inBam=bamFile, outBed=rmDupedBedFile))
      	     add(new Pyicoclip(inBed = rmDupedBedFile, outBed = pyicoclipResults, regions = genicRegionsLocation(genome) ))

	     var ripseeker = new RipSeeker
	     ripseeker.inBam = bamFile
	     ripseeker.outBed = swapExt(bamFile, ".bam", ".ripseeker.bed")
     
	     var piranha = new Piranha
	     piranha.inBam = bamFile
	     piranha.outBed = swapExt(bamFile, ".bam", ".piranha.bed")
	     
	     var clipClassic = new ClipClassic
	     clipClassic.inBam = bamFile
	     clipClassic.species = genome
	     clipClassic.out_file = kasey_output
	     add(clipClassic)
	     add(new ClipAnalysis(bamFile, kasey_output, genome, kasey_output_metrics, AS_Structure = asStructureLocation(genome), 
	     			  genome_location = genomeLocation(genome), phastcons_location = phastconsLocation(genome), 
	     			  gff_db = gffDbLocation(genome)))
	     

	     var miclip = new MiClip
	     miclip.inBam = bamFile
	     miclip.genome = genomeLocation(genome)
	     miclip.outBed = swapExt(bamFile, ".bam", ".miclip.bed")
	     add(miclip)

	     var pipeclip = new PipeClip
	     pipeclip.inBam = bamFile
	     pipeclip.species = genome
	     pipeclip.outBed = swapExt(bamFile, ".bam", ".pipeclip.bed")
	     add(pipeclip)

	     add(ripseeker, piranha)

      	     //add(new IDR(inBam = bamFile, species = genome, genome = chromSizeLocation(genome), outResult = IDRResult, premRNA = premRNA))

      	     add(new countTags(input = bamFile, index = bamIndex, output = countFile, a = exonLocation(genome)))
      	     add(new singleRPKM(input = countFile, output = RPKMFile, s = genome))



  }

  def script() {

    val fileList = QScriptUtils.createArgsFromFile(input)
    var posTrackHubFiles: List[File] = List()
    var negTrackHubFiles: List[File] = List()
   
    for ((groupName, valueList) <- (fileList groupBy (_._3))) {
    	 var combinedBams : Seq[File] = List()
	 var genome: String = valueList(0)._2 
	 for (item : Tuple3[File, String, String] <- valueList) {
	     var fastq_file: File = item._1

      	     var replicate: String = item._3

      	     val noPolyAFastq = swapExt(fastq_file, ".fastq", ".polyATrim.fastq")
      	     val noPolyAReport = swapExt(noPolyAFastq, ".fastq", ".metrics")

      	     val noAdapterFastq = swapExt(noPolyAFastq, ".fastq", ".adapterTrim.fastq")
      	     val adapterReport = swapExt(noAdapterFastq, ".fastq", ".metrics")

      	     val filteredFastq = swapExt(noAdapterFastq, ".fastq", ".rmRep.fastq")
      	     val filterd_results = swapExt(filteredFastq, ".fastq", ".metrics")

      	     val samFile = swapExt(filteredFastq, ".fastq", ".sam")
	     val rgBamFile = swapExt(samFile, ".sam", ".rg.bam")
      	     var sortedrgBamFile = swapExt(rgBamFile, ".bam", ".sorted.bam")

      	     val NRFFile = swapExt(sortedrgBamFile, ".bam", ".NRF.metrics")

      	     val rmDupedBamFile = swapExt(sortedrgBamFile, ".bam", ".rmDup.bam")
      	     val rmDupedMetricsFile = swapExt(rmDupedBamFile, ".bam", ".metrics")

      	     val sortedrmDupedBamFile = swapExt(rmDupedBamFile, ".bam", ".sorted.bam")
      	     val indexedBamFile = swapExt(sortedrmDupedBamFile, "", ".bai")
	     
	     if(!fromBam) { 
      	     	add(new FastQC(inFastq = fastq_file))
		//filters out adapter reads
	      	add(cutadapt(fastq_file = fastq_file, noAdapterFastq = noAdapterFastq, adapterReport = adapterReport, adapter = adapter ) )
	     
		add(filterRepetitiveRegions(noAdapterFastq, filterd_results, filteredFastq))
      	     	add(new FastQC(filteredFastq))
      	     	add(star(input = filteredFastq, output = samFile, genome_location = starGenomeLocation(genome)))
      	     	add(addOrReplaceReadGroups(samFile, rgBamFile))
	     	add(sortSam(rgBamFile, sortedrgBamFile, SortOrder.coordinate))

      	     } else {
	        sortedrgBamFile = fastq_file
	     }
      	     add(new CalculateNRF(inBam = sortedrgBamFile, genomeSize = chromSizeLocation(genome), outNRF = NRFFile))
      	     add(new RemoveDuplicates(sortedrgBamFile, rmDupedBamFile, rmDupedMetricsFile))
      	     add(sortSam(rmDupedBamFile, sortedrmDupedBamFile, SortOrder.coordinate))
      	     add(new samtoolsIndexFunction(sortedrmDupedBamFile, indexedBamFile))
	     combinedBams = combinedBams ++ List(sortedrmDupedBamFile)
	     downstream_analysis(sortedrmDupedBamFile, indexedBamFile, genome)
    	     }

	//Only run combination if we are combining more than one thing...
	if(combinedBams.size > 1 && groupName != "null") {
		var mergedBams = new File(groupName + ".bam")
		val mergedIndex = swapExt(mergedBams, "", ".bai")
		add(new samtoolsIndexFunction(mergedBams, mergedIndex))
		add(samtoolsMergeFunction(combinedBams, mergedBams))
		downstream_analysis(mergedBams, mergedIndex, genome)
	}
        }
    }
}
