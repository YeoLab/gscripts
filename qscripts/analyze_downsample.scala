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

    this.nCoresRequest = Option(8)
    this.wallTime = Option((24 * 60 * 60).toLong)
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


def call_peaks(bam: File, input_bam: File, genome: String) {

  val bam_peaks = swapExt(bam, ".bam", ".peaks.bed")
  val bam_peaks_norm = swapExt(bam_peaks, ".bed", ".norm.bed")
  
  val bam_peaks_compressed = swapExt(bam_peaks_norm, ".bed", ".compressed.bed")
  val bam_peaks_fixed = swapExt(bam_peaks_compressed, ".bed", ".l2Fixed.bed")

  add(new clipper(in = bam, genome = genome, out = bam_peaks, isPremRNA = premRNA, reverse=reverse_strand))
  add(new InputNorm(bam_peaks, bam, input_bam, bam_peaks_norm, bam_peaks_compressed, bam_peaks_fixed))
}

def split_and_call_peaks(bam: File, input_bam: File, genome: String) {
  val bam_01 = swapExt(bam, ".bam", ".01.bam")
  val bam_02 = swapExt(bam, ".bam", ".02.bam")
  val bam_03 = swapExt(bam, ".bam", ".03.bam")
  val bam_04 = swapExt(bam, ".bam", ".04.bam")
  val bam_05 = swapExt(bam, ".bam", ".05.bam")
  val bam_06 = swapExt(bam, ".bam", ".06.bam")
  val bam_07 = swapExt(bam, ".bam", ".07.bam")
  val bam_08 = swapExt(bam, ".bam", ".08.bam")
  val bam_09 = swapExt(bam, ".bam", ".09.bam")
  
  add(new DownsampleBam(bam, bam_01, bam_02, bam_03, bam_04, bam_05, bam_06, bam_07, bam_08, bam_09))
  
  call_peaks(bam_01, input_bam, genome)
  call_peaks(bam_02, input_bam, genome)
  call_peaks(bam_03, input_bam, genome)
  call_peaks(bam_04, input_bam, genome)
  call_peaks(bam_05, input_bam, genome)
  call_peaks(bam_06, input_bam, genome)
  call_peaks(bam_07, input_bam, genome)
  call_peaks(bam_08, input_bam, genome)
  call_peaks(bam_09, input_bam, genome)

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
      split_and_call_peaks(rep1, input, genome)
      split_and_call_peaks(rep2, input, genome)

    }
  }
}






