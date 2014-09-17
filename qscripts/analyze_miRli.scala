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
import scala.util.parsing.combinator._
import scala.io.Source


class Analyze_mirli_CLIPSeq extends QScript {

    // Script args
    @Input(doc = "input file")
    var input: File = _

    @Argument(doc = "bam_counts")
    var bam_counts: String = "bam_counts.txt"

    var bam_counts_file = new File(bam_counts)

    @Argument(doc = "adapter to trim")
    var adapter: List[String] = Nil

    @Argument(doc = "RBP binding modality is premRNA")
    var premRNA: Boolean = false

    @Argument(doc = "read ids have randomers")
    var barcoded: Boolean = false

    @Argument(doc = "use alpha version of STAR")
    var alpha: Boolean = false

    @Argument(doc = "miRbase file", required=true)
    var mirbase: String = "none"

    @Argument(doc = "base quality miniumum", required=false)
    var qual_cut: java.lang.Integer  = 13

    @Argument(doc = "read length miniumum", required=false)
    var len_cut: java.lang.Integer  = 36

    @Argument(doc = "genome", required=true)
    var genome: String = "ce10"

    @Argument(doc = "regionsToMask", required=false)
    var regionsToMask_str: String = "/home/jbrought/scratch/mirpipe/background_removal_test/20140828.WS240.background_RNA.galaxy1.bed"

    var regionsToMask: File = new File(regionsToMask_str)

    //https://gist.github.com/paradigmatic/3437345
    object FASTA {

      case class Entry( description: String, sequence: String )

      def fromFile( fn: String ): List[Entry] = {
        val lines = io.Source.fromFile(fn).getLines.mkString("\n")
        fromString( lines )
      }

      def fromString( input: String ): List[Entry] =
        Parser.parse(input)

      private object Parser extends RegexParsers {

        lazy val header = """>.*""".r ^^ { _.tail.trim }
        lazy val seqLine = """[^>].*""".r ^^ { _.trim }
        lazy val sequence = rep1( seqLine ) ^^ { _.mkString }
        lazy val entry = header ~ sequence ^^ {
          case h ~ s => Entry(h,s)
        }

        lazy val entries = rep1( entry )
        def parse( input: String ): List[Entry]  = {
          parseAll( entries, input ) match {
            case Success( es , _ ) => es
            case x: NoSuccess =>  throw new Exception(x.toString)
          }
        }
      }
    }
    // end https://gist.github.com/paradigmatic/3437345

    class FastqGroom(@Input inFastq: File, @Output outFastq: File) extends CommandLineFunction {
        override def shortDescription = "FastqGroom"
        def commandLine = "fastq_quality_trimmer " +
        required("-t", qual_cut) +
        required("-l", len_cut) +
        required("-i", inFastq) +
        required("-o", outFastq)
    }

    class SplitBy(@Input inFastq: File,
                  @Output splitFastq: File,
                  @Argument miRNA_seq: String,
                  @Argument miRNA_id: String ) extends CommandLineFunction {

        override def shortDescription = "SplitBy"
        def commandLine = "miR_splitter.py " +
                          required("--miRNA_seq", miRNA_seq) +
                          required("--miRNA_id", miRNA_id) +
                          required("--fastq", inFastq) +
                          required("--output", splitFastq)
    }

    case class clipper(in: File, out: File, genome: String, isPremRNA: Boolean ) extends Clipper
    {
        this.inBam = in
        this.outBed = out
        this.species = genome
        this.premRNA = isPremRNA
        this.superlocal = superlocal
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
        def commandLine =   "bamToBed " +
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


    case class maskRegions(@Input maskMe: File, maskWith: File, finished: File, @Output maskedOut: File) extends CommandLineFunction {

        override def shortDescription = "maskRepeats"
        def commandLine = "bedtools window -v -header -w 50" +
        required("-a", maskMe) +
        required("-b", maskWith) +
        + " > " + maskedOut

    }

    case class countBamToGenes(@Input bamFiles: List[File], geneBedFile: File, @Output bamCounts: File) extends CommandLineFunction{

        override def shortDescription = "countBamToGenes"
        def commandLine = "bedtools multicov -bams " +
        repeat(bamFiles) +
        " -s -D -bed " + geneBedFile " >> " + bamCounts
    }

    case class mergeBam(@Input bamFile: File, @Output mergedBed: File) extends CommandLineFunction{

        override def short Description = "mergeBam"
        def commandLine = "bamToBed -i " + bamFile + " -splitD > " + mergedBed

    }
    def downstream_analysis(bamFile : File, bamIndex: File, genome : String) = {

        val bedGraphFilePos = swapExt(bamFile, ".bam", ".pos.bg")
        val bedGraphFilePosNorm = swapExt(bedGraphFilePos, ".pos.bg", ".norm.pos.bg")
        val bigWigFilePos = swapExt(bedGraphFilePosNorm, ".bg", ".bw")

        val bedGraphFileNeg = swapExt(bamFile, ".bam", ".neg.bg")
        val bedGraphFileNegNorm = swapExt(bedGraphFileNeg, ".neg.bg", ".norm.neg.bg")
        val bigWigFileNeg = swapExt(bedGraphFileNeg, ".bg", ".normal.bw")
        val bedGraphFileNegInverted = swapExt(bedGraphFileNegNorm, "neg.bg", "neg.t.bg")
        val bigWigFileNegInverted = swapExt(bedGraphFileNegInverted, ".t.bg", ".bw")

        val clipper_output = swapExt(bamFile, ".bam", ".peaks.bed")
        val fixed_clipper_output = swapExt(clipper_output, ".bed", ".fixed.bed")
        val bigBed_output = swapExt(fixed_clipper_output, ".bed", ".bb")
        val clipper_output_metrics = swapExt(clipper_output, ".bed", ".metrics")

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
        add(new BedGraphToBigWig(bedGraphFileNegNorm, chromSizeLocation(genome), bigWigFileNeg))
        add(new NegBedGraph(inBedGraph = bedGraphFileNegNorm, outBedGraph = bedGraphFileNegInverted))
        add(new BedGraphToBigWig(bedGraphFileNegInverted, chromSizeLocation(genome), bigWigFileNegInverted))

        add(new clipper(in = bamFile, genome = genome, out = clipper_output, isPremRNA = premRNA))

        add(new FixScores(inBed = clipper_output, outBed = fixed_clipper_output))

        add(new BedToBigBed(inBed = fixed_clipper_output, genomeSize = chromSizeLocation(genome), outBigBed = bigBed_output))

        add(new ClipAnalysis(bamFile, clipper_output, genome, clipper_output_metrics,
        regions_location = regionsLocation(genome), AS_Structure = asStructureLocation(genome),
        genome_location = genomeLocation(genome), phastcons_location = phastconsLocation(genome),
        gff_db = gffDbLocation(genome), bw_pos=bigWigFilePos, bw_neg=bigWigFileNeg))

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

        //add(ripseeker, piranha, clipClassic)
        //add(new IDR(inBam = bamFile, species = genome, genome = chromSizeLocation(genome), outResult = IDRResult, premRNA = premRNA))
        //add(new countTags(input = bamFile, index = bamIndex, output = countFile, a = exonLocation(genome)))
        //add(new singleRPKM(input = countFile, output = RPKMFile, s = genome))
    }

    def script() = {

        val mirbase_entries = FASTA.fromFile( mirbase )

        var fastq_files = QScriptUtils.createSeqFromFile(input)
        var finished_bams_files: List[File] = List()



        for (fastq_filename <- fastq_files) {

            var fastq_file = new File(fastq_filename)
            val noPolyAFastq = swapExt(fastq_file, ".fastq", ".polyATrim.fastq")
            val noPolyAReport = swapExt(noPolyAFastq, ".fastq", ".metrics")
            val noAdapterFastq = swapExt(noPolyAFastq, ".fastq", ".adapterTrim.fastq")
            val adapterReport = swapExt(noAdapterFastq, ".fastq", ".metrics")

            add(cutadapt(fastq_file = fastq_file, noAdapterFastq = noAdapterFastq,
                         adapterReport = adapterReport, adapter = adapter))

            val groomedFastq: File = swapExt(noAdapterFastq, ".fastq", ".groom.fastq")

            add(new FastqGroom(noAdapterFastq, groomedFastq))
            add(new FastQC(inFastq = groomedFastq))

            for( mir <- mirbase_entries ) {

                var mir_id = mir.description
                var mir_seq = mir.sequence

                val miRFastq = swapExt(groomedFastq, ".fastq", "." + mir_id + ".fastq")

                add(new SplitBy(groomedFastq, miRFastq, mir_seq, mir_id))
                add(new FastQC(miRFastq))

                val samFile = swapExt(miRFastq, ".fastq", ".sam")
                add(star(input = miRFastq, output = samFile, genome_location = starGenomeLocation(genome)))

                val rgBamFile = swapExt(samFile, ".sam", ".rg.bam")
                add(addOrReplaceReadGroups(samFile, rgBamFile))

                val sortedrgBamFile = swapExt(rgBamFile, ".bam", ".sorted.bam")
                add(sortSam(rgBamFile, sortedrgBamFile, SortOrder.coordinate))

                val NRFFile = swapExt(sortedrgBamFile, ".bam", ".NRF.metrics")
                add(new CalculateNRF(inBam = sortedrgBamFile, genomeSize = chromSizeLocation(genome), outNRF = NRFFile))

                val rmDupedBamFile = swapExt(sortedrgBamFile, ".bam", ".rmDup.bam")
                val rmDupedMetricsFile = swapExt(rmDupedBamFile, ".bam", ".metrics")
                add(new RemoveDuplicates(sortedrgBamFile, rmDupedBamFile, rmDupedMetricsFile))

                val sortedrmDupedBamFile = swapExt(rmDupedBamFile, ".bam", ".sorted.bam")
                add(sortSam(rmDupedBamFile, sortedrmDupedBamFile, SortOrder.coordinate))

                val indexedBamFile = swapExt(sortedrmDupedBamFile, "", ".bai")
                add(new samtoolsIndexFunction(sortedrmDupedBamFile, indexedBamFile))

                val maskedSortedrmDupedBamFile = swapExt(sortedrmDupedBamFile, ".bam", ".masked")
                add(new maskRegions(sortedrmDupedBamFile, regionsToMask, finished_bams_file,
                                    maskedSortedrmDupedBamFile))

                finished_bams_files = finished_bams_files ++ List(maskedSortedrmDupedBamFile)

                val indexedMaskedBamFile = swapExt(maskedSortedrmDupedBamFile, "", ".bai")
                add(new samtoolsIndexFunction(maskedSortedrmDupedBamFile, indexedMaskedBamFile))

                val mergedMaskedBamFile = swapExt(maskedSortedrmDupedBamFile, "", ".merged.bed")
                add(new mergeBam(maskedSortedrmDupedBamFile, mergedMaskedBamFile))
                val mergedMaskedBamFileMetrics = swapExt(mergedMaskedBamFile, ".bed", ".metrics")


                add(new ClipAnalysis(maskedSortedrmDupedBamFile, mergedMaskedBamFile, genome, mergedMaskedBamFileMetrics,
                    regions_location = regionsLocation(genome), AS_Structure = asStructureLocation(genome),
                    genome_location = genomeLocation(genome), phastcons_location = phastconsLocation(genome),
                    gff_db = gffDbLocation(genome), bw_pos=bigWigFilePos, bw_neg=bigWigFileNeg))



                //downstream_analysis(sortedrmDupedBamFile, indexedBamFile, genome)
            }
        }

        countBamToGenes(finished_bams_file, genicRegionsLocation(genome), bam_counts_file)
    }
}