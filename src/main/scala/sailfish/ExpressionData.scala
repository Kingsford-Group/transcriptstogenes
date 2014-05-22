import scala.collection.mutable.ArrayBuffer
import io.Source
import scalax.io._
import processing.Processor
import scalax.file.Path

class ExpressionDatum(val name: String, val length: Int, var quants: Array[Double]) { //}val tpm: Double, val rpkm: Double, val kpkm: Double, val nreads: Double) {
  override def toString() = { 
    //s"${name}\t${length}\t${tpm}\t${rpkm}\t${kpkm}\t${nreads}"    
    s"""${name}\t${length}\t${quants.mkString("\t")}"""
  }
}

object ExpressionData {
  
  /**
  * Read a sailfish format expression file and return a tuple containing the 
  * header information and the expression dataset
  */
  def fromFile(fname: String): (String, ExpressionData) = {
    
    var dset:ArrayBuffer[ExpressionDatum] = new ArrayBuffer
    var header = ""
    val fullName = fname.replace("~", System.getProperty("user.home"))
    val path: Path = Path.fromString(fullName).toAbsolute
    val lineIt = io.Source.fromFile(path.path).getLines
    // Parse each line as a header or expression line
    lineIt.foreach({ l =>
      if (!l.startsWith("#")) { 
        val toks = l.split("\t")
        dset += new ExpressionDatum(toks(0), toks(1).toDouble.toInt, toks.tail.tail.map{ _.toDouble } )//.toDouble, toks(3).toDouble, toks(4).toDouble, toks(5).toDouble)
      } else { 
        header += l
        header += "\n"
      }
    })
    
    println(s"Expression file contained ${dset.size} lines")
    // Return the concatenated header string and an ExpressionData instance containing the parsed
    // expression data
    (header, new ExpressionData(dset.toArray))
   }
  
  /**
  * Write the header ("header") and expression information ("expd") to a file with the name
  * given in "outfile".
  */
  def toFile(outFile: String, header: String, expd: ExpressionData) { 
    val fullName = outFile.replace("~", System.getProperty("user.home"))
    val path: Path = Path(fullName).toAbsolute
    path.createFile(failIfExists=false)
    
    val output = Resource.fromFile(path.path)
    var ctr = 0
    for{
       processor <- output.outputProcessor
       out = processor.asOutput
    }{ 
      out.write(header)
      expd.dset.foreach({ datum => ctr += 1; out.write(datum + "\n") })
    }
    println(s"wrote out ${ctr} records\n")

  }
}

/**
* Represents and expression dataset.  Its a simple row-oriented format, which is just a collection of 
* ExpressionDatum instances.
*/
class ExpressionData(val dset:Array[ExpressionDatum]) { 
}

