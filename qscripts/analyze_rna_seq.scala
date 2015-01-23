package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.picard.{ ReorderSam, SortSam, AddOrReplaceReadGroups, MarkDuplicates }
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.util.YeowareUtils._
import org.broadinstitute.sting.queue.extensions.yeo._

class AnalyzeRNASeq extends QScript {
  // Script argunment
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

  @Argument(doc = "start processing run from bam file")
  var fromBam: Boolean = false

  @Argument(doc = "location to place trackhub (must have the rest of the track hub made before starting script)")
  var location: String = "rna_seq"

 @Argument(doc="reads are single ended", shortName = "single_end", fullName = "single_end", required = false)
 var singleEnd: Boolean = true
    
  case class sortSam(inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam {
    override def shortDescription = "sortSam"
    this.input :+= inSam
    this.output = outBam
    this.sortOrder = sortOrderP
    this.createIndex = true
  }

  case class filterRepetitiveRegions(noAdapterFastq: File, filteredResults: File, filteredFastq: File) extends FilterRepetitiveRegions {
       override def shortDescription = "FilterRepetitiveRegions"

       this.inFastq = noAdapterFastq
       this.outCounts = filteredResults
       this.outNoRep = filteredFastq
       this.isIntermediate = true
  }

  case class genomeCoverageBed(input: File, outBed : File, cur_strand: String, species: String) extends GenomeCoverageBed {
        this.inBam = input
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
    this.tags_annotation = exonLocation(species)
    this.flip = flipped
  }
  
  case class star(input: File, output: File, stranded : Boolean, paired : File = null, species: String) extends STAR {
    this.inFastq = input
    
    if(paired != null) {
      this.inFastqPair = paired
    }
    
    this.outSam = output
    //intron motif should be used if the data is not stranded
    this.intronMotif = stranded
    this.genome = starGenomeLocation(species) 
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
  
  case class cutadapt(fastqFile: File, noAdapterFastq: File, adapterReport: File, adapter: List[String]) extends Cutadapt{
    override def shortDescription = "cutadapt"
    
    this.inFastq = fastqFile
    this.outFastq = noAdapterFastq
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
  

  
  def stringentJobs(fastqFile: File) : File = {
    
    // run if stringent
      val noPolyAFastq = swapExt(fastqFile, ".fastq", ".polyATrim.fastq")
    val noPolyAReport = swapExt(noPolyAFastq, ".fastq", ".metrics")
    val noAdapterFastq = swapExt(noPolyAFastq, ".fastq", ".adapterTrim.fastq")
    val adapterReport = swapExt(noAdapterFastq, ".fastq", ".metrics")
    val filteredFastq = swapExt(noAdapterFastq, ".fastq", ".rmRep.fastq")
    val filtered_results = swapExt(filteredFastq, ".fastq", ".metrics")
    //filters out adapter reads
    add(cutadapt(fastqFile = fastqFile, noAdapterFastq = noAdapterFastq, 
		 adapterReport = adapterReport, 
		 adapter = adapter))
    
    
    add(filterRepetitiveRegions(noAdapterFastq, filtered_results, filteredFastq))
    add(new FastQC(filteredFastq))
    
    return filteredFastq
  }
  
  def makeBigWig(inBam: File, species: String): (File, File) = {
    
    val bedGraphFilePos = swapExt(inBam, ".bam", ".pos.bg")
    val bedGraphFilePosNorm = swapExt(bedGraphFilePos, "pos.bg", ".norm.pos.bg")
    val bigWigFilePosNorm = swapExt(bedGraphFilePosNorm, ".bg", ".bw")
    
    val bedGraphFileNeg = swapExt(inBam, ".bam", ".neg.bg")
    val bedGraphFileNegNorm = swapExt(bedGraphFileNeg, "neg.bg", ".norm.neg.bg")
    val bedGraphFileNegNormInverted = swapExt(bedGraphFileNegNorm, ".bg", ".t.bg")
    val bigWigFileNegNormInverted = swapExt(bedGraphFileNegNormInverted, ".t.bg", ".bw")
    
    add(new genomeCoverageBed(input = inBam, outBed = bedGraphFilePos, cur_strand = "+", species = species))
    add(new NormalizeBedGraph(inBedGraph = bedGraphFilePos, inBam = inBam, outBedGraph = bedGraphFilePosNorm))
    add(new BedGraphToBigWig(bedGraphFilePosNorm, chromSizeLocation(species), bigWigFilePosNorm))
    
    add(new genomeCoverageBed(input = inBam, outBed = bedGraphFileNeg, cur_strand = "-", species = species))
    add(new NormalizeBedGraph(inBedGraph = bedGraphFileNeg, inBam = inBam, outBedGraph = bedGraphFileNegNorm))
    add(new NegBedGraph(inBedGraph = bedGraphFileNegNorm, outBedGraph = bedGraphFileNegNormInverted))
    add(new BedGraphToBigWig(bedGraphFileNegNormInverted, chromSizeLocation(species), bigWigFileNegNormInverted))
    return (bigWigFileNegNormInverted, bigWigFilePosNorm)
    
  }
  

  case class samtoolsMergeFunction(inBams: Seq[File], outBam: File) extends SamtoolsMergeFunction {
    override def shortDescription = "samtoolsMerge"
    this.inputBams = inBams
    this.outputBam = outBam
  }

  def downstream_analysis(bamFile: File, bamIndex: File, singleEnd: Boolean, species: String) : File = {
    val NRFFile = swapExt(bamFile, ".bam", ".NRF.metrics")      
    val countFile = swapExt(bamFile, "bam", "count")
    val RPKMFile = swapExt(countFile, "count", "rpkm")
    val oldSpliceOut = swapExt(bamFile, "bam", "splices")
    val misoOut = swapExt(bamFile, "bam", "miso")
    val rnaEditingOut = swapExt(bamFile, "bam", "editing")
    
    add(new CalculateNRF(inBam = bamFile, genomeSize = chromSizeLocation(species), outNRF = NRFFile))
    
    val (bigWigFilePos: File, bigWigFileNeg: File) = makeBigWig(bamFile, species = species)
        
    add(new countTags(input = bamFile, index = bamIndex, output = countFile, species = species))	
    add(new singleRPKM(input = countFile, output = RPKMFile, s = species))
    
    add(oldSplice(input = bamFile, out = oldSpliceOut, species = species))
    add(new Miso(inBam = bamFile, indexFile = bamIndex, species = species, pairedEnd = false, output = misoOut))
    add(new RnaEditing(inBam = bamFile, snpEffDb = species, snpDb = snpDbLocation(species), genome = genomeLocation(species), flipped=flipped, output = rnaEditingOut))
    return oldSpliceOut
  }

  def script() {
    
    
    
    val fileList = QScriptUtils.createArgsFromFile(input)
    var posTrackHubFiles: List[File] = List()
    var negTrackHubFiles: List[File] = List()
    var splicesFiles: List[File] = List()
    var bamFiles: List[File] = List()
    var speciesList: List[String] = List()
    for ((groupName, valueList) <- (fileList groupBy (_._3))) {
      var combinedBams : Seq[File] = List()
      var genome: String = valueList(0)._2 

      for (item : Tuple3[File, String, String] <- valueList) {
	var fastqFiles = item._1.toString().split(""";""")
	var species = item._2
	var fastqFile: File = new File(fastqFiles(0))
	var fastqPair: File = null
	var singleEnd = true
	var samFile: File = null
	if(!fromBam) {
      	  if (fastqFiles.length == 2){
            singleEnd = false
            fastqPair = new File(fastqFiles(1))
	    add(new FastQC(inFastq = fastqPair))
      	  }
          
     	  add(new FastQC(inFastq = fastqFile))
	  
      	  var filteredFastq: File = null
      	  if(strict && fastqPair == null) {
       	    filteredFastq = stringentJobs(fastqFile)
      	  } else {
            filteredFastq = fastqFile
      	  }
	  samFile = swapExt(filteredFastq, ".fastq", ".sam")
      	 
    sailfishOutDir = swapExt(swapExt(fastqFile, ".gz", ""), ".fastq", ".sailfish")
    sailfishShScript = swapExt(sailfishOutDir, ".sailfish", ".sailfish.sh")
    sailfishIndex = SailfishGenomeIndexLocation(species)

	  if(fastqPair != null) {
	        //if paired
            add(new sailfish(inFastq=filteredFastq, index=sailfishIndex,
              outDir=sailfishOutDir, stranded=!not_stranded, shScript=sailfishShScript, 
              inFastqPair=fastqPair)

            add(new star(filteredFastq, samFile, not_stranded, fastqPair, species = species))
      	  } else { //unpaired
            add(new sailfish(inFastq=filteredFastq, index=sailfishIndex,
              outDir=sailfishOutDir, stranded=!not_stranded, shScript=sailfishShScript)

            add(new star(filteredFastq, samFile, not_stranded, species = species))

	  }
	  
	  // run regardless of stringency
      	  
	} else {
      	  samFile = new File(fastqFiles(0))
	}
	
	val sortedBamFile = swapExt(samFile, ".sam", ".sorted.bam")
	val rgSortedBamFile = swapExt(sortedBamFile, ".bam", ".rg.bam")
	val indexedBamFile = swapExt(rgSortedBamFile, "", ".bai")
	
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






