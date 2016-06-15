\name{LibraryReportDDT}
\alias{LibraryReportDDT}
\alias{LibraryReportKelthane}
\title{Mass Spectral Library Reports}
\description{Generate reports for the two mass spectral libraries: DDT technical mixture and Kelthane-Dicofol technical mixture.}
\usage{LibraryReportDDT(spectra = SpecData,
                     metadata = MetaDataDDT,
                     pdfFile = "SpecLibDDT2015.pdf",
                     pdfTitle = "SpecLibDDT2015 Library",
                     xMin = 40)
}

\arguments{
   \item{spectra}{data frame containing the mass spectra.}
    \item{metadata}{data frame containing the metadata.}
    \item{pdfFile}{character string. The file name of the PDF report.}
    \item{pdfTitle}{character string. Title string to embed in the PDF file.} 
    \item{xMin}{numeric value. The minimum m/z in the plot.}
}
\details{See the OrgMassSpecR - Mass Spectral Libraries vignette for an
  overview and examples.
}
\value{A PDF report of the mass spectra.  If PDF = FALSE, the report is
  generated using the graphics device.
}

