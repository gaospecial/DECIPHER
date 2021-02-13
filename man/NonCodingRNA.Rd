\name{NonCodingRNA}
\alias{NonCodingRNA_Archaea}
\alias{NonCodingRNA_Bacteria}
\alias{NonCodingRNA_Eukarya}
\docType{data}
\title{
NonCoding Models for Common Non-Coding RNA Families
}
\description{
Pre-trained with \code{NonCoding} models for common RNA families found in genomes from organisms belonging to each domain of life.
}
\usage{data("NonCodingRNA_Archaea")}
\details{
A set of \code{NonCoding} models contained in a named list.  Models were built from up to 1000 representative sequences per non-coding RNA family.
}
\source{
Models were built from sequences belonging to families in tRNADB-CE (\url{http://trna.ie.niigata-u.ac.jp/cgi-bin/trnadb/index.cgi}) or Rfam (\url{http://rfam.xfam.org}).
}
\examples{
data(NonCodingRNA_Archaea)
data(NonCodingRNA_Bacteria)
data(NonCodingRNA_Eukarya)
names(NonCodingRNA_Bacteria)
head(NonCodingRNA_Bacteria)
}
\keyword{datasets}
