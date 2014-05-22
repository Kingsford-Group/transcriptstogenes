import io.Source
import scala.collection.mutable.{ArrayBuffer, OpenHashMap => HMap, HashSet => HSet}
import org.rogach.scallop._
import scalax.io._
import processing.Processor
import scalax.file.Path

object StratifyBySimilarity {
  
  /**
  * Configuration options
  */
  class Conf(arguments: Seq[String]) extends ScallopConf(arguments) { 
      version("0.1.0")
      banner("Transcript-level to gene-level expression conversion")
      footer("========================================")
      val gtfFile = opt[String](required=true, descr="GTF file containing the gene-to-transcript mapping")
      //val expFile = opt[String](required=true, descr="Transcript-level expression file")
      val blastFile = opt[String](required=true, descr="File containing transcript-level BLAST")
      val resFile = opt[String](required=false, descr="File where gene-level expression should be written")
  }
      
  /**
  * Parses a line from a GTF file and returns the "transcript_id" and "gene_id" attributes
  */
  def getAttributes(l:String) = {
    val tab = """\t""".r
    val space = """\s+""".r
    val attributes = l.split("\t").last.split(";") //ab.split(l).last.split(";")
    var gid = ""
    var tid = ""
    val amap = attributes.foreach({ a => 
        val Array(k,v) = a.trim.split(" ")//space.split(a.trim)
        k match {
           case "gene_id" => gid = v.stripPrefix("\"").stripSuffix("\"")
           case "transcript_id" => tid = v.stripPrefix("\"").stripSuffix("\"")
           case _ => 
        }
    })
    (gid, tid)
  }

  object BLASTEntry { 
    def fromTokens(toks: Array[String]) = { 
        require(toks.length == 12)
        new BLASTEntry(toks(0), toks(1), toks(2).toFloat, toks(3).toInt, toks(4).toInt,
                       toks(5).toInt, toks(6).toInt, toks(7).toInt, toks(8).toInt, toks(9).toInt,
                       toks(10).toFloat, toks(11).toFloat)
    }  
  }
  
  class BLASTEntry(val query: String, val subject: String, val percentID: Float, val alnLen: Int, val mismatches: Int,
                   val gapOpenings: Int, val queryStart: Int, val queryEnd: Int, val subjectStart: Int, val subjectEnd: Int,
                   val eValue: Float, val bitScore: Float) { 
  
  
        def hitToString():String = { 
          s"${query}, %id=${percentID}, alnLength=${alnLen}, e-value=${eValue}, bit-score=${bitScore}" 
        }
  }
  
  def parseBLASTFile(blastPath: Path, tgm: HMap[String, String]) = { 
    val lineIt = io.Source.fromFile(blastPath.path).getLines
    val geneMap: HMap[String, ArrayBuffer[BLASTEntry]] = new HMap
    var i = 0
    lineIt.foreach({ l => 
      val toks = l.split("\t")
      val queryGene = tgm(toks(0)); toks(0) = queryGene
      val subjectGene = tgm(toks(1)); toks(1) = subjectGene
      val entry = BLASTEntry.fromTokens(toks)
      
      if (!geneMap.contains(subjectGene)) { 
        geneMap(subjectGene) = new ArrayBuffer[BLASTEntry]
      }       
      geneMap(subjectGene) += entry
      if (i % 10000 == 0) { print(s"Processed ${i} BLAST entries\r\r") }
      i += 1
    })
    
    geneMap.par.foreach({ case (key: String, value: ArrayBuffer[BLASTEntry]) =>
      val seenSet:HSet[String] = new HSet 
      geneMap(key) = value.filter( { x => 
          val keep = (x.query != key) && !seenSet.contains(x.query); seenSet += x.query; keep 
      }).sortBy( _.eValue )
    }) 
    
    println
    geneMap
  }
  
  /**
  * Main driver function.  This program aggregates transcript-level expression estimates
  * output by Sailfish into gene-level expression estimates.
  */
  def main(args: Array[String]) = {
    
    val conf = new Conf(args)
    
    // Transcript => Gene and Gene => Transcript maps
    val tgm:HMap[String, String] = new HMap
    val gtm:HMap[String, HSet[String]] = new HMap
    
    // Parse the GTF file
    val fn = conf.gtfFile().replace("~", System.getProperty("user.home"))
    val path: Path = Path.fromString(fn).toAbsolute
    val lineIt = io.Source.fromFile(path.path).getLines
    var i = 0
    // For each line, add the gene <=> transcript association
    lineIt.foreach({ l => 
      val (gid, tid) = getAttributes(l) 
      tgm(tid) = gid
      gtm(gid) = gtm.get(gid) getOrElse new HSet[String]()
      gtm(gid) += tid
      i += 1
      if (i % 100000 == 0) { print(s"\r\rline ${i}") }
    })
    println()
  
    println(s"There were ${tgm.size} transcripts and ${gtm.size} genes")
    println("Parsing BLAST results")
     
    val blastFn = conf.blastFile().replace("~", System.getProperty("user.home"))
    val blastPath: Path = Path.fromString(blastFn).toAbsolute
    val geneSim = parseBLASTFile(blastPath, tgm)
   
    val outFile = Resource.fromFile(conf.resFile().replace("~", System.getProperty("user.home")))
   
    val observedGenes:HSet[String] = new HSet
    for { 
        outProcessor <- outFile.outputProcessor
        out = outProcessor.asOutput
      } { 
        geneSim.foreach({ case (key:String, value:ArrayBuffer[BLASTEntry]) =>
          if (value.length > 0) { 
             observedGenes.add(key)
             value.foreach{ v =>
               out.write(s"${key} => ${v.hitToString}\n") 
             }
          }
        })
        
        val remainingGenes = gtm.keySet &~ observedGenes
        remainingGenes.foreach({ g =>
             out.write(s"${g} => NONE, %id=0.0, alnLength=0, e-value=100.0, bit-score=0.0, size(0)\n")
        })
        
      }
     
    /*
    // Parse the sailfish file, returning the header 
    // (information about the quantification run) and the expression dataset.
    val (header, ed) = ExpressionData.fromFile(conf.expFile())
    
    // For each element of the expression dataset, map it from transcript to 
    // gene land.  We keep the length of the gene as the length of the longest
    // transcript.  The expression values simply get summed up.
    val expmap:HMap[String, ExpressionDatum] = new HMap
    ed.dset.foreach({ datum => 
      val gid = tgm(datum.name)
      val v = expmap.get(gid) getOrElse new ExpressionDatum("", 0, 0.0, 0.0, 0.0, 0.0)
      expmap(gid) = new ExpressionDatum(gid, math.max(v.length, datum.length), v.tpm + datum.tpm, v.rpkm + datum.rpkm, v.kpkm + datum.kpkm, v.nreads + datum.nreads)
    })
    
    // Create an expression dataset from the gene-level expression values
    val geneExpression = new ExpressionData(expmap.values.toArray)
    
    println(s"There were ${geneExpression.dset.size} elements at the gene level.")
    
    // Write the gene-level expression dataset to file
    print(s"Writing gene-level results to ${conf.resFile()} . . . ")
    ExpressionData.toFile(conf.resFile(), header, geneExpression)
    println("done")
    */
  }
}
