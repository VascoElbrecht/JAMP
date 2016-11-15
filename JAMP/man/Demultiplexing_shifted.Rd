\name{Demultiplexing_shifted}
\alias{Demultiplexing_shifted}
\title{
Demultiplexing PE illumina HTS data based on inline barcodes}
\description{
R script to demultiplex samples tagged from both sides using the inline tagging (Frame shift fusion primers, see Elbrecht & Leese 2015).  
}
\usage{
}
\arguments{
  \item{file1}{First illumina reads. Can be one or several fastq files with "gz" and "bz2" compression or directly as "fastq".}
  \item{file2}{Seqcond illumina reads. See \code{file1}.}
  \item{tags}{Information on the tags used (in csv format). Common fusion primers have already tables bundled with JAMP. Just write the respective primer names; \code{"BF_BR"} (Elbrecht et al. 2016).}
  \item{combinations}{Provide information on used tagging combinations and respecitve sample names as a csv table. See Tutorial folder of the package fo an example.}
  \item{md5}{Calculates md5 checksums for the raw data and displays the information in the R interface and \code{log.txt}.}
}
\details{
nope}
\value{
}
\references{
Elbrecht & Leese (2015). Can DNA-Based Ecosystem Assessments Quantify Species Abundance? Testing Primer Bias and Biomass—Sequence Relationships with an Innovative Metabarcoding Protocol. PlosONE
Elbrecht, Peinert, Leese (2016). Validation and development of freshwater invertebrate metabarcoding COI primers for Environmental Impact Assessment. PeerJ PrePrints
}
\author{
}
\note{
}


\seealso{
\code{\link{NA}}
}
\examples{}