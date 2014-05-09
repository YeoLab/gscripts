package org.broadinstitute.sting.queue.qscripts

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import net.sf.samtools.SAMFileHeader.SortOrder
import org.broadinstitute.sting.utils.exceptions.UserException
import org.broadinstitute.sting.commandline.Hidden
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.yeo._
import org.broadinstitute.sting.queue.util.TsccUtils._
import scala.io.Source

class DoStuff extends QScript {

  // runs a list of command-line commands using GATK Queue
  // --input is a file containing commands that "work" on the command line
  // --ncores is the number of cores required for each job
  // --jobname is the name of the submitted jobs

  @Input(doc = "input file")
  var input: File = _

  @Argument(doc = "n cores for jobs", required=false)
  var ncores: Int = 1

  @Argument(doc = "job name", required=false)
  var jobname: String = "eval"

  case class evalNth(@Input fileIn: File, @Argument lineN: String) extends CommandLineFunction{
       override def shortDescription = jobname + "." + lineN
       this.nCoresRequest = Option(ncores)
       def commandLine = "eval `head -n " + lineN + " " + fileIn + " | tail -n 1`"
  }

  def looper(@Input input: File){
    for(item <- 1 to Source.fromFile(input).getLines.size){
    	add(new evalNth(fileIn=input, lineN=item.toString()))
    }
  }

  def script() {
    looper(input)
  }

}