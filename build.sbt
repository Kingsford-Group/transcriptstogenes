import AssemblyKeys._

assemblySettings

jarName in assembly := "TranscriptsToGenes.jar"

mainClass in assembly := Some("TranscriptToGene")

name :="TranscriptToGene"

scalaVersion :="2.10.3"

version :="1.0"

libraryDependencies += "com.github.scala-incubator.io" %% "scala-io-core" % "0.4.2"

libraryDependencies += "com.github.scala-incubator.io" %% "scala-io-file" % "0.4.2"

libraryDependencies += "org.rogach" %% "scallop" % "0.9.4"
