\name{nifti.pbcor}
\alias{nifti.pbcor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Parcel-based correlation between two NIfTI objects}
\description{Calculates the correlation between two NIfTI objects (e.g., brain maps) across random parcels.}
\usage{
nifti.pbcor(nifti1, nifti2, mask_coords, n.parcels = 90, n.parcellations = 50,
            kmeans.iter.max = 30, kmeans.nstart = 1,
            kmeans.algorithm = "Hartigan-Wong", kmeans.trace = FALSE,
            cor.alternative = "two.sided", cor.method = "pearson",
            cor.exact = NULL, cor.conf.level = 0.95, cor.continuity = FALSE,
            verbose = TRUE)
}
\arguments{
  \item{nifti1}{an object of class \code{"nifti"} with 3D (e.g., the first brain map).}
  \item{nifti2}{an object of class \code{"nifti"} with 3D (e.g., the second brain map).}
  \item{mask_coords}{a mask matrix prepared by \link{nifti.pbcor_mask}.}
  \item{n.parcels}{the number of parcels to divide the mask into. The default is 90, approximately the minimum number of resels in the meta-analytic maps of gray matter differences between individuals with mental disorders and healthy controls in Fortea et al. (see below).}
  \item{n.parcellations}{the number of parcellations, i.e., how many times the mask will be parcellated.}
  \item{kmeans.iter.max}{(optional) argument \code{iter.max} passed to \code{\link{kmeans}}}
  \item{kmeans.nstart}{(optional) argument \code{nstart} passed to \code{\link{kmeans}}}
  \item{kmeans.algorithm}{(optional) argument \code{algorithm} passed to \code{\link{kmeans}}}
  \item{kmeans.trace}{(optional) argument \code{trace} passed to \code{\link{kmeans}}}
  \item{cor.alternative}{(optional) argument \code{alternative} passed to \code{\link{cor.test}}}
  \item{cor.method}{(optional) argument \code{method} passed to \code{\link{cor.test}}}
  \item{cor.exact}{(optional) argument \code{exact} passed to \code{\link{cor.test}}}
  \item{cor.conf.level}{(optional) argument \code{conf.level} passed to \code{\link{cor.test}}}
  \item{cor.continuity}{(optional) argument \code{continuity} passed to \code{\link{cor.test}}}
  \item{verbose}{(optional) logical, whether to print some messages during execution.}
}
\details{This approach resolves the problems of voxel-based correlations, where contiguous voxels are non-independent, by randomly segmenting the mask (e.g., the brain) into parcels and calculating the Pearson correlation coefficient across these parcels. Additionally, it circumvents the limitations of atlas-based correlations, where estimates depend on the specific atlas used, by repeating the random parcellation and correlation multiple times and selecting the median estimate.}
\value{The parcel-based correlation between the two images across parcellations. The attribute "parcellations.cor.test" has the \code{\link{cor.test}} results obtained in the different random parcellations.}
\references{Fortea et al., under review.}
\author{Joaquim Radua}
\seealso{\code{\link{nifti.pbcor_mask}}, \code{\link[oro.nifti]{readNIfTI}}}
\examples{
library(oro.nifti)

# Path of the example files (with large voxels to ensure the example runs quickly)
mask_path = system.file("examples", "mask.nii.gz", package = "nifti.pbcor")
img1_path = system.file("examples", "img1.nii.gz", package = "nifti.pbcor")
img2_path = system.file("examples", "img2.nii.gz", package = "nifti.pbcor")

# Prepare the mask
mask = nifti.pbcor_mask(readNIfTI(mask_path))

# Conduct the parcel-based correlation
nifti.pbcor(readNIfTI(img1_path), readNIfTI(img2_path), mask)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
