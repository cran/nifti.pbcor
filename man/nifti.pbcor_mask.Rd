\name{nifti.pbcor_mask}
\alias{nifti.pbcor_mask}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prepare the mask for a parcel-based correlation between NIfTI images}
\description{This function creates the mask required to conduct parcel-based correlations with \code{\link{nifti.pbcor}}.}
\usage{
nifti.pbcor_mask(nifti, verbose = TRUE)
}
\arguments{
  \item{nifti}{an object of class \code{"nifti"} with 3D (e.g., the brain mask).}
  \item{verbose}{(optional) logical, whether to print some messages during execution.}
}
\details{This function converts a NIfTI mask into the appropriate format for \code{\link{nifti.pbcor}}.}
\value{A matrix with the coordinates of the mask}
\author{Joaquim Radua}
\seealso{\code{\link{nifti.pbcor}}, \code{\link[oro.nifti]{readNIfTI}}}
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
