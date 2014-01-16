package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.picard.{ ReorderSam, SortSam, AddOrReplaceReadGroups, MarkDuplicates }
import org.broadinstitute.sting.queue.extensions.samtools._
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.util.TsccUtils._
import org.broadinstitute.sting.queue.extensions.yeo._

class AnalyzeRNASeq extends QScript {
  // Script argunment
  @Input(doc = "input file or txt file of input files")
  var input: File = _

  @Argument(doc = "species (hg19...)")
  var species: String = _

  @Argument(doc = "adapter to trim")
  var adapter: List[String] = Nil

  @Argument(doc = "flipped", required = false)
  var flipped: String = _

  @Argument(doc = "not stranded", required = false)
  var not_stranded: Boolean = false

  @Argument(doc = "strict triming run")
  var strict: Boolean = false

  @Argument(doc = "location to place trackhub (must have the rest of the track hub made before starting script)")
  var location: String = "rna_seq"
    
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

  case class genomeCoverageBed(input: File, outBed : File, cur_strand: String) extends GenomeCoverageBed {
        this.inBam = input
        this.genomeSize = chromSizeLocation(species)
        this.bedGraph = outBed
        this.strand = cur_strand
        this.split = true
  }

  
  case class oldSplice(input: File, out : File) extends OldSplice {
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

  case class countTags(input: File, index: File, output: File) extends CountTags {
	this.inBam = input
	this.outCount = output
	this.tags_annotation = exonLocation(species)
	this.flip = flipped
  }

  case class star(input: File, output: File, stranded : Boolean, paired : File = null) extends STAR {
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

  case class parseOldSplice(inSplice : List[File]) extends ParseOldSplice {
       override def shortDescription = "ParseOldSplice"
       this.inBam = inSplice
       this.in_species = species
  }
 
  case class samtoolsIndexFunction(input: File, output: File) extends SamtoolsIndexFunction {
       override def shortDescription = "indexBam"
       this.bamFile = input
       this.bamFileIndex = output
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


  case class makeRNASeQC(input: List[File], output: File) extends MakeRNASeQC {
       this.inBam = input
       this.out = output

}
  case class runRNASeQC(in : File, out : String, single_end : Boolean) extends RunRNASeQC { 
       this.input = in
       this.gc = gcLocation(species)
       this.output = out
       this.genome = genomeLocation(species)
       this.gtf = gffLocation(species)
       this.singleEnd = single_end
}
  

 @Argument(doc="reads are single ended", shortName = "single_end", fullName = "single_end", required = false)
 var singleEnd: Boolean = true

def stringentJobs(fastq_file: File) : File = {

        // run if stringent
      val noPolyAFastq = swapExt(fastq_file, ".fastq", ".polyATrim.fastq")
      val noPolyAReport = swapExt(noPolyAFastq, ".fastq", ".metrics")
      val noAdapterFastq = swapExt(noPolyAFastq, ".fastq", ".adapterTrim.fastq")
      val adapterReport = swapExt(noAdapterFastq, ".fastq", ".metrics")
      val filteredFastq = swapExt(noAdapterFastq, ".fastq", ".rmRep.fastq")
      val filtered_results = swapExt(filteredFastq, ".fastq", ".metrics")
            //filters out adapter reads
      add(cutadapt(fastq_file = fastq_file, noAdapterFastq = noAdapterFastq, 
          adapterReport = adapterReport, 
          adapter = adapter))
          
          
      add(filterRepetitiveRegions(noAdapterFastq, filtered_results, filteredFastq))
      add(new FastQC(filteredFastq))

      return filteredFastq
}
def script() {

    
    
    val fileList = QScriptUtils.createArgsFromFile(input)
    var trackHubFiles: List[File] = List()
    var splicesFiles: List[File] = List()
    var bamFiles: List[File] = List()
    
    for (item : Tuple3[File, String, String] <- fileList) {
      var fastq_file: File = item._1
      var fastqPair: File = null
      var singleEnd = true
      if (item._2 != "null"){
        singleEnd = false
        fastqPair = new File(item._2)
	add(new FastQC(inFastq = fastqPair))
      }
      
      add(new FastQC(inFastq = fastq_file))

      var filteredFastq: File = null
      if(strict && item._2 == "null") {
        filteredFastq = stringentJobs(fastq_file)
      } else {
        filteredFastq = fastq_file
      }


      // run regardless of stringency
      val samFile = swapExt(filteredFastq, ".fastq", ".sam")
      val sortedBamFile = swapExt(samFile, ".sam", ".sorted.bam")
      val rgSortedBamFile = swapExt(sortedBamFile, ".bam", ".rg.bam")
      val indexedBamFile = swapExt(rgSortedBamFile, "", ".bai")
      
      val NRFFile = swapExt(rgSortedBamFile, ".bam", ".NRF.metrics")

      val bedGraphFilePos = swapExt(rgSortedBamFile, ".bam", ".pos.bg")
      val bedGraphFilePosNorm = swapExt(bedGraphFilePos, ".pos.bg", ".norm.pos.bg")
      val bigWigFilePos = swapExt(bedGraphFilePosNorm, ".bg", ".bw")

      val bedGraphFileNeg = swapExt(rgSortedBamFile, ".bam", ".neg.bg")
      val bedGraphFileNegNorm = swapExt(bedGraphFileNeg, "neg.bg", ".norm.neg.bg")
      val bedGraphFileNegInverted = swapExt(bedGraphFileNegNorm, ".bg", ".t.bg")
      val bigWigFileNeg = swapExt(bedGraphFileNegInverted, ".t.bg", ".bw")

      val countFile = swapExt(rgSortedBamFile, "bam", "count")
      val RPKMFile = swapExt(countFile, "count", "rpkm")

      val oldSpliceOut = swapExt(rgSortedBamFile, "bam", "splices")
      
      //add bw files to list for printing out later

      trackHubFiles = trackHubFiles ++ List(bigWigFileNeg, bigWigFilePos)
      splicesFiles = splicesFiles ++ List(oldSpliceOut)      
      bamFiles = bamFiles ++ List(rgSortedBamFile)

      if(item._2 != "null") { //if paired	
        add(new star(filteredFastq, samFile, not_stranded, fastqPair))
      } else { //unpaired
        add(new star(filteredFastq, samFile, not_stranded))
      }
      
      add(new sortSam(samFile, sortedBamFile, SortOrder.coordinate))
      add(addOrReplaceReadGroups(sortedBamFile, rgSortedBamFile))
      add(new samtoolsIndexFunction(rgSortedBamFile, indexedBamFile))
      add(new CalculateNRF(inBam = rgSortedBamFile, genomeSize = chromSizeLocation(species), outNRF = NRFFile))
      
      add(new genomeCoverageBed(input = rgSortedBamFile, outBed = bedGraphFilePos, cur_strand = "+"))
      add(new NormalizeBedGraph(inBedGraph = bedGraphFilePos, inBam = rgSortedBamFile, outBedGraph = bedGraphFilePosNorm))
      add(new BedGraphToBigWig(bedGraphFilePosNorm, chromSizeLocation(species), bigWigFilePos))

      add(new genomeCoverageBed(input = rgSortedBamFile, outBed = bedGraphFileNeg, cur_strand = "-"))
      add(new NormalizeBedGraph(inBedGraph = bedGraphFileNeg, inBam = rgSortedBamFile, outBedGraph = bedGraphFileNegNorm))
      add(new NegBedGraph(inBedGraph = bedGraphFileNegNorm, outBedGraph = bedGraphFileNegInverted))
      add(new BedGraphToBigWig(bedGraphFileNegInverted, chromSizeLocation(species), bigWigFileNeg))
	
      add(new countTags(input = rgSortedBamFile, index = indexedBamFile, output = countFile))	
      add(new singleRPKM(input = countFile, output = RPKMFile, s = species))

      add(oldSplice(input = rgSortedBamFile, out = oldSpliceOut))
    }
    add(new MakeTrackHub(trackHubFiles, location, species))
    add(parseOldSplice(splicesFiles))
    var rnaseqc = new File("rnaseqc.txt")
    add(new makeRNASeQC(input = bamFiles, output = rnaseqc))
    add(new runRNASeQC(in = rnaseqc, out = "rnaseqc", single_end = singleEnd))
  }
}






