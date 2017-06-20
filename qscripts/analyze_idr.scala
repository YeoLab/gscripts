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
     
  @Argument(doc = "RBP binding modality is premRNA")
  var premRNA: Boolean = false

  @Argument(doc = "reads are on opposite strand")
  var reverse_strand: Boolean = false

  case class samtoolsIndexFunction(input: File, output: File) extends SamtoolsIndexFunction {
       override def shortDescription = "indexBam"
       this.bamFile = input
       this.bamFileIndex = output
  }


  case class clipper(in: File, out: File, genome: String, isPremRNA: Boolean, reverse: Boolean) extends Clipper
  {

    this.wallTime = Option((48 * 60 * 60).toLong)
    this.inBam = in
    this.outBed = out
    this.species = genome
    this.premRNA = isPremRNA
    this.superlocal = superlocal
    this.reverse_strand = reverse
  }
  

  case class samtoolsMergeFunction(inBams: Seq[File], outBam: File) extends SamtoolsMergeFunction {
    override def shortDescription = "samtoolsMerge"
    this.inputBams = inBams
    this.outputBam = outBam
  }

 /* case class splitBam(inBam: File, out01: File, out02: File) extends SplitBam {
    override def shortDescription = "splitBam"
    this.bam = inBam
    this.bam01 = out01
    this.bam02 = out02
  }

*/


def split_and_call_peaks(bam: File, input_bam: File, genome: String): (File, File) = {
  val bam_01 = swapExt(bam, ".bam", ".01.bam")
  val bam_02 = swapExt(bam, ".bam", ".02.bam")
  add(new SplitBam(bam, bam_01, bam_02))
  
  //val bam_01_index = swapExt(bam_01, "", ".bai")
  //val bam_02_index = swapExt(bam_02, "", ".bai")
  //add(new samtoolsIndexFunction(bam_01, bam_01_index))
  //add(new samtoolsIndexFunction(bam_02, bam_02_index))

  val bam_01_peaks = swapExt(bam_01, ".bam", ".peaks.bed")
  val bam_02_peaks = swapExt(bam_02, ".bam", ".peaks.bed")

  val bam_01_peaks_norm = swapExt(bam_01_peaks, ".bed", ".norm.bed")
  val bam_02_peaks_norm = swapExt(bam_02_peaks, ".bed", ".norm.bed")

  val bam_01_peaks_compressed = swapExt(bam_01_peaks_norm, ".bed", ".compressed.bed")
  val bam_02_peaks_compressed = swapExt(bam_02_peaks_norm, ".bed", ".compressed.bed")

  val bam_01_peaks_fixed = swapExt(bam_01_peaks_compressed, ".bed", ".l2Fixed.bed")
  val bam_02_peaks_fixed = swapExt(bam_02_peaks_compressed, ".bed", ".l2Fixed.bed")

  add(new clipper(in = bam_01, genome = genome, out = bam_01_peaks, isPremRNA = premRNA, reverse=reverse_strand))
  add(new InputNorm(bam_01_peaks, bam_01, input_bam, bam_01_peaks_norm, bam_01_peaks_compressed, bam_01_peaks_fixed))
  
  add(new clipper(in = bam_02, genome = genome, out = bam_02_peaks, isPremRNA = premRNA, reverse=reverse_strand))
  add(new InputNorm(bam_02_peaks, bam_02, input_bam, bam_02_peaks_norm, bam_02_peaks_compressed, bam_02_peaks_fixed))
  
  return (bam_01_peaks_fixed, bam_02_peaks_fixed)

}
  def script() {
  
    val fileList = QScriptUtils.createArgsFromFile(input)
  
    for (item : Tuple7[File, String, String, String, String, String, String] <- fileList) {  
      var rep1 = item._1.toString()
      var rep2 = item._2
      var input = item._3
      var genome = item._4
      var groupName = item._5
    
      //split rep1 and rep2
      var (rep1_01, rep1_02) = split_and_call_peaks(rep1, input, genome)
      var (rep2_01, rep2_02) =  split_and_call_peaks(rep2, input, genome)
      var rep1_idr = new IDR_V2(rep1_01, rep1_02, new File(groupName + ".rep1.txt"))
      rep1_idr.wallTime = Option((24 * 60 * 60).toLong)
      add(rep1_idr)
      var rep2_idr = new IDR_V2(rep2_01, rep2_02, new File(groupName + ".rep2.txt"))
      rep2_idr.wallTime = Option((24 * 60 * 60).toLong)
      add(rep2_idr)
      //merge and run IDR
      var mergedBams = new File(groupName + ".merged.bam")
      val combinedBams = List(rep1, rep2) 
      add(samtoolsMergeFunction(combinedBams, mergedBams))
      var (combined_01, combined_02) = split_and_call_peaks(mergedBams, input, genome)
      var combined_idr = new IDR_V2(combined_01, combined_02, new File(groupName + ".combined.txt"))
      combined_idr.wallTime = Option((24 * 60 * 60).toLong)
      add(combined_idr)
      
      // Run IDR on the real peaks this is super hacky because I'll need to check peak results after, make sure they are exactly the same as eric's results
      var rep1_real = swapExt(rep1.getParent, rep1, ".bam", ".peaks.bed")
      var rep2_real = swapExt(rep2.getParent, rep2, ".bam", ".peaks.bed")
      
      val bam_01_peaks_norm = swapExt(rep1_real, ".bed", ".norm.bed")
      val bam_02_peaks_norm = swapExt(rep2_real, ".bed", ".norm.bed")

      val bam_01_peaks_compressed = swapExt(bam_01_peaks_norm, ".bed", ".compressed.bed")
      val bam_02_peaks_compressed = swapExt(bam_02_peaks_norm, ".bed", ".compressed.bed")

      val bam_01_peaks_fixed = swapExt(bam_01_peaks_compressed, ".bed", ".l2Fixed.bed")
      val bam_02_peaks_fixed = swapExt(bam_02_peaks_compressed, ".bed", ".l2Fixed.bed")

      add(new InputNorm(rep1_real, rep1, input, bam_01_peaks_norm, bam_01_peaks_compressed, bam_01_peaks_fixed))
      add(new InputNorm(rep2_real, rep2, input, bam_02_peaks_norm, bam_02_peaks_compressed, bam_02_peaks_fixed))

      var normal_idr = new IDR_V2(bam_01_peaks_fixed, bam_02_peaks_fixed, new File(groupName + ".txt"))
      normal_idr.wallTime = Option((24 * 60 * 60).toLong)
      add(normal_idr)
    }
  }
}






