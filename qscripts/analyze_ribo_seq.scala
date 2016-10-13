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
class AnalyzeRNASeq extends QScript {
  // Script argument
  @Input(doc = "input file or txt file of input files")
  var input: File = _

  @Argument(doc = "adapter to trim")
  var adapter: List[String] = Nil

  @Argument(doc = "flipped", required = false)
  var flipped: String = "none"

  @Argument(doc = "not stranded", required = false)
  var not_stranded: Boolean = false

  @Argument(doc = "strict triming run")
  var strict: Boolean = false

  @Argument(doc = "start processing from bam file")
  var fromBam: Boolean = false

  @Argument(doc = "location to place trackhub (must have the rest of the track hub made before starting script)")
  var location: String = "rna_seq"

  case class sortSam(inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam {
    override def shortDescription = "sortSam"
    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP
    this.createIndex = true
  }

  
  case class genomeCoverageBed(input: File, outBed : File, cur_strand: String, species: String) extends GenomeCoverageBed {
        this.inBed = input
        this.genomeSize = chromSizeLocation(species)
        this.bedGraph = outBed
        this.strand = cur_strand
        this.split = true
  }

  case class oldSplice(input: File, out : File, species: String) extends OldSplice {
        this.inBam = input
        this.out_file = out
        this.in_species = species
        this.splice_type = List("MXE", "SE")
        this.flip = (flipped == "flip")
  }
  case class singleRPKM(input: File, output: File, s: String) extends SingleRPKM {
	this.inCount = input
	this.outRPKM = output
  }

  case class countTags(input: File, index: File, output: File, species: String) extends CountTags {
	this.inBam = input
	this.outCount = output
	this.tags_annotation =  exonLocation(species)
	this.flip = flipped
  }


  case class star(input: File,                                                                                                                                                                                                                                                                                                     
                  output: File,                                                                                                                                                                                                                                                                                                    
                  stranded: Boolean,                                                                                                                                                                                                                                                                                               
                  genome_location: String,                                                                                                                                                                                                                                                                                         
                  paired: File = null,
                  fastq_out: File = null
                  ) extends STAR {
  
    this.inFastq = input
  
    if (paired != null) {
      this.inFastqPair = paired
    }
    
    this.outSam = output
    //intron motif should be used if the data is not stranded                                                                                                                                                                                                                                                                      
    this.intronMotif = stranded
    this.genome = genome_location
    this.outFastq = fastq_out
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
       this.RGID = "foo"

  }

  case class parseOldSplice(inSplice : List[File], species: String) extends ParseOldSplice {
       override def shortDescription = "ParseOldSplice"
       this.inBam = inSplice
       this.in_species = species
  }
 
  case class samtoolsIndexFunction(input: File, output: File) extends SamtoolsIndexFunction {
       override def shortDescription = "indexBam"
       this.bamFile = input
       this.bamFileIndex = output
  }
  case class cutadapt(fastqFile: File, noAdapterFastq: File, adapterReport: File, adapter: List[String], pairedFile: File = null, pairedOut: File = null, a_adapters: List[String] = Nil, A_adapters: List[String] = Nil, g_adapters: List[String] = Nil, G_adapters: List[String] = Nil) extends Cutadapt {
    override def shortDescription = "cutadapt"
  
    this.inFastq = fastqFile
    this.outFastq = noAdapterFastq

    if(pairedFile != null) {
      this.inPair = pairedFile
      this.outPair = pairedOut
    }
    this.three_prime = a_adapters
    this.three_prime2 = A_adapters
    this.five_prime = g_adapters
    this.five_prime2 = G_adapters
    this.report = adapterReport
    this.anywhere = adapter ++ List("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
      "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
    this.overlap = 5
    this.length = 18
    this.quality_cutoff = 6
    this.isIntermediate = true
  } 



  case class makeRNASeQC(input: List[File], output: File) extends MakeRNASeQC {
       this.inBam = input
       this.out = output

}
  case class runRNASeQC(in : File, out : String, single_end : Boolean, species: String) extends RunRNASeQC {
       this.input = in
       this.gc = gcLocation(species)
       this.output = out
       this.genome = genomeLocation(species)
       this.gtf = gffLocation(species)
       this.singleEnd = single_end
}

  case class fastQC(fastq: File, dir: String= null) extends FastQC {
    this.inFastq = fastq
    this.outDir = dir
  } 

  def stringentJobs(fastqFile: File, pairedFile: File = null): (File, File) = {

    // run if stringent                                                                                                                                                                                                                      
    val noGzFastq = swapExt(fastqFile, ".gz", "")
    val noAdapterFastq = swapExt(noGzFastq, ".fastq", ".polyATrim.adapterTrim.fastq")
    val outRep = swapExt(noAdapterFastq, ".fastq", ".rep.bam")

    val filteredFastq = swapExt(outRep, "", "Unmapped.out.mate1")
    val adapterReport = swapExt(noAdapterFastq, ".fastq", ".metrics")

    var noGzPair: File = null
    var noAdapterPair: File = null
    var filteredPair: File = null

    if (pairedFile != null) {
     noGzPair = swapExt(pairedFile, ".gz", "")
     noAdapterPair = swapExt(noGzPair, ".fastq", ".polyATrim.adapterTrim.fastq")
     filteredPair = swapExt(outRep, "", "Unmapped.out.mate2")
    }

    //filters out adapter reads                                                                                                                                                                                                             
    add(cutadapt(fastqFile = fastqFile,
                 noAdapterFastq = noAdapterFastq,
                 adapterReport = adapterReport,
		 adapter = adapter,
             
                 pairedFile = pairedFile,
                 pairedOut = noAdapterPair))
    
    add(star(input = noAdapterFastq,
             paired = noAdapterPair,
             output = outRep,
             stranded = not_stranded,
             genome_location = "/projects/ps-yeolab/genomes/RepBase18.05.fasta/STAR_fixed",
             fastq_out = filteredFastq
	     ))

    var countRepetitiveRegions = new CountRepetitiveRegions
    countRepetitiveRegions.inBam = outRep
    countRepetitiveRegions.outFile = swapExt(outRep, ".rep.bam", ".rmRep.metrics")
    add(countRepetitiveRegions)
    
    add(new fastQC(filteredFastq, dir=qSettings.runDirectory))

    if (filteredPair != null) {
      add(new fastQC(filteredPair, dir=qSettings.runDirectory))
    }

    return (filteredFastq, filteredPair)
  }
  
  case class samtoolsMergeFunction(inBams: Seq[File], outBam: File) extends SamtoolsMergeFunction {
    override def shortDescription = "samtoolsMerge"
    this.inputBams = inBams
    this.outputBam = outBam
  }
  
  def makeBigWig(inBed: File, inBam: File, species: String): (File, File) = {

    val bedGraphFilePos = swapExt(inBed, ".bed", ".pos.bg")
    val bedGraphFilePosNorm = swapExt(bedGraphFilePos, ".bg", ".norm.bg")
    val bigWigFilePosNorm = swapExt(bedGraphFilePosNorm, ".bg", ".bw")
    val bigWigFilePos = swapExt(bedGraphFilePos, ".bg", ".bw")
    
    val bedGraphFileNeg = swapExt(inBed, ".bed", ".neg.bg")
    val bedGraphFileNegInverted = swapExt(bedGraphFileNeg, ".bg", ".t.bg")
    val bedGraphFileNegInvertedNorm = swapExt(bedGraphFileNegInverted, ".bg", ".norm.bg")
    val bigWigFileNegNorm = swapExt(bedGraphFileNegInvertedNorm, ".bg", ".bw")
    val bigWigFileNeg = swapExt(bedGraphFileNegInverted, ".bg", ".bw")
    
    add(new genomeCoverageBed(input = inBed, outBed = bedGraphFilePos, cur_strand = "+", species = species))
    add(new NormalizeBedGraph(inBedGraph = bedGraphFilePos, inBam = inBam, outBedGraph = bedGraphFilePosNorm))
    add(new BedGraphToBigWig(bedGraphFilePosNorm, chromSizeLocation(species), bigWigFilePosNorm))
    add(new BedGraphToBigWig(bedGraphFilePos, chromSizeLocation(species), bigWigFilePos))
    
    add(new genomeCoverageBed(input = inBed, outBed = bedGraphFileNeg, cur_strand = "-", species = species))
    add(new NegBedGraph(inBedGraph = bedGraphFileNeg, outBedGraph = bedGraphFileNegInverted))
    add(new NormalizeBedGraph(inBedGraph = bedGraphFileNegInverted, inBam = inBam, outBedGraph = bedGraphFileNegInvertedNorm))
    add(new BedGraphToBigWig(bedGraphFileNegInvertedNorm, chromSizeLocation(species), bigWigFileNegNorm))
    add(new BedGraphToBigWig(bedGraphFileNegInverted, chromSizeLocation(species), bigWigFileNeg))
    return (bigWigFileNegNorm, bigWigFilePosNorm)

  }

  def downstream_analysis(bamFile: File, bamIndex: File, singleEnd: Boolean, species: String) : File = {
    val bedFileAdjusted = swapExt(bamFile, ".bam", ".adj.bed")
    val NRFFile = swapExt(bamFile, ".bam", ".NRF.metrics")
    
    val countFile = swapExt(bamFile, "bam", "count")
    val RPKMFile = swapExt(bamFile, "count", "rpkm")
    
    val oldSpliceOut = swapExt(bamFile, "bam", "splices")
    val misoOut = swapExt(bamFile, "bam", "miso")

    add(new CalculateNRF(inBam = bamFile, genomeSize = chromSizeLocation(species), outNRF = NRFFile))
    
    add(new RiboSeqCoverage(inBam = bamFile, outBed = bedFileAdjusted))
    val (bigWigFilePos: File, bigWigFileNeg: File) = makeBigWig(bedFileAdjusted, bamFile, species = species)
    //posTrackHubFiles = posTrackHubFiles ++ List(bigWigFilePos)
    //negTrackHubFiles = negTrackHubFiles ++ List(bigWigFileNeg)
    
    add(new countTags(input = bamFile, index = bamIndex, output = countFile, species = species))
    add(new singleRPKM(input = countFile, output = RPKMFile, s = species))
    
    add(oldSplice(input = bamFile, out = oldSpliceOut, species = species))
    add(new Miso(inBam = bamFile, indexFile= bamIndex, species = species, pairedEnd = !singleEnd, output = misoOut))

    return oldSpliceOut
  }

  def script() {
    val fileList = QScriptUtils.createArgsFromFile(input)
    var posTrackHubFiles: List[File] = List()
    var negTrackHubFiles: List[File] = List()
    var splicesFiles: List[File] = List()
    var speciesList: List[File] = List()
    var bamFiles: List[File] = List()
    var singleEnd = true
    for ((groupName, valueList) <- (fileList groupBy (_._3))) {
      var combinedBams : Seq[File] = List()
      var genome: String = valueList(0)._2

      for (item : Tuple7[File, String, String, String, String, String, String] <- valueList) {
	
	var fastq_files = item._1.toString().split(""";""")
	var species = item._2
	var replicate: String = item._3
	var samFile: File = null
	if(!fromBam) {
      	  var fastq_file: File = new File(fastq_files(0))
      	  var fastqPair: File = null
      	  var singleEnd = true
	  
     	  if (fastq_files.length == 2){
            singleEnd = false
            fastqPair = new File(fastq_files(1))
            add(new fastQC(fastqPair, qSettings.runDirectory))
      	  }
	  
	  add(new fastQC(fastq_file, qSettings.runDirectory))
	  
     
	  val (filteredFastq, filteredPair) = stringentJobs(fastq_file, fastqPair)
     	  
	  
	
	  // run regardless of stringency
      	  samFile = swapExt(filteredFastq, ".fastq", ".sam")
	
	  add(star(input = filteredFastq,
		   paired = fastqPair,
		   output = samFile,
		   stranded = not_stranded,
		   genome_location = starGenomeLocation(species)
		 ))
       	  
	  
	} else {
          samFile = new File(fastq_files(0))
	}
	
	val sortedBamFile = swapExt(samFile, ".sam", ".sorted.bam")
	val rgSortedBamFile = swapExt(sortedBamFile, ".bam", ".rg.bam")
	val indexedBamFile = swapExt(rgSortedBamFile, "", ".bai")      
	//add bw files to list for printing out later
	
	
	bamFiles = bamFiles ++ List(rgSortedBamFile)
	speciesList = speciesList ++ List(species)
      
	add(new sortSam(samFile, sortedBamFile, SortOrder.coordinate))
	add(addOrReplaceReadGroups(sortedBamFile, rgSortedBamFile))
	add(new samtoolsIndexFunction(rgSortedBamFile, indexedBamFile))
	combinedBams = combinedBams ++ List(rgSortedBamFile)
	var oldSpliceOut = downstream_analysis(rgSortedBamFile, indexedBamFile, singleEnd, genome)
	splicesFiles = splicesFiles ++ List(oldSpliceOut)      
      }  
      
      if(groupName != "null") {
	var mergedBams = new File(groupName + ".bam")
	val mergedIndex = swapExt(mergedBams, "", ".bai")
	add(new samtoolsIndexFunction(mergedBams, mergedIndex))
	add(samtoolsMergeFunction(combinedBams, mergedBams))
	var oldSpliceOut = downstream_analysis(mergedBams, mergedIndex, singleEnd, genome)
	splicesFiles = splicesFiles ++ List(oldSpliceOut)   
      }
      
    }
  
    def tuple2ToList[T](t: (T,T)): List[T] = List(t._1, t._2)
    for ((species, files) <- posTrackHubFiles zip negTrackHubFiles zip speciesList groupBy {_._2}) {
      var trackHubFiles = files map {_._1} map tuple2ToList reduceLeft {(a,b) => a ++ b}
      add(new MakeTrackHub(trackHubFiles, location=location + "_" + species, genome=species))
    }
    
    for ((species, files) <- speciesList zip splicesFiles groupBy {_._1}) {
      add(parseOldSplice(files map {_._2}, species = species))
    }
  
    for ((species, files) <- speciesList zip bamFiles groupBy {_._1}) {
      var rnaseqc = new File("rnaseqc_" + species + ".txt")
      add(new makeRNASeQC(input = files map {_._2}, output = rnaseqc))
      add(new runRNASeQC(in = rnaseqc, out = "rnaseqc_" + species, single_end = singleEnd, species = species))
    }
  }
}






