import io.Source
import scala.collection.mutable.{OpenHashMap => HMap, HashSet => HSet}
import org.rogach.scallop._
import scalax.file.Path

object TranscriptToGene{
  
  /**
  * Configuration options
  */
  class Conf(arguments: Seq[String]) extends ScallopConf(arguments) { 
      version("0.1.0")
      banner("Transcript-level to gene-level expression conversion")
      footer("========================================")
      val gtfFile = opt[String](required=true, descr="GTF file containing the gene-to-transcript mapping")
      val expFile = opt[String](required=true, descr="Transcript-level expression file")
      val resFile = opt[String](required=true, descr="File where gene-level expression should be written")
      val minVal = opt[Double](required=false, descr="Minimum expression value to consider non-zero", default=Option(0.01))
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
    
    // Parse the sailfish file, returning the header 
    // (information about the quantification run) and the expression dataset.
    val (header, ed) = ExpressionData.fromFile(conf.expFile())
    
    // For each element of the expression dataset, map it from transcript to 
    // gene land.  We keep the length of the gene as the length of the longest
    // transcript.  The expression values simply get summed up.
    val expmap:HMap[String, ExpressionDatum] = new HMap
    ed.dset.foreach({ datum =>
      if (tgm.contains(datum.name)) {
          val gid = tgm(datum.name)
          val v = expmap.get(gid) getOrElse new ExpressionDatum("", 0, Array.fill(datum.quants.size){ 0.0 })
          v.quants = v.quants.map({ q => if (q < conf.minVal()) { 0.0 } else { q } })
          expmap(gid) = new ExpressionDatum(gid, math.max(v.length, datum.length), datum.quants.zip(v.quants).map( x => x._1 + x._2 ) )
       } else {
          println(s"GTF did not contain entry for ${datum.name}; expression was ${datum.quants(0)}")
      }
    })
   
    val numQuant = ed.dset(0).quants.size
    gtm.foreach({ case (key: String, value: HSet[String]) =>
      if (!expmap.contains(key)) { 
          expmap(key) = new ExpressionDatum(key, 0, Array.fill(numQuant){ 0.0 }) 
      }  
    })
    // Create an expression dataset from the gene-level expression values
    val geneExpression = new ExpressionData(expmap.values.toArray)
    
    println(s"There were ${geneExpression.dset.size} elements at the gene level.")
    
    // Write the gene-level expression dataset to file
    print(s"Writing gene-level results to ${conf.resFile()} . . . ")
    ExpressionData.toFile(conf.resFile(), header, geneExpression)
    println("done")
  }
}
