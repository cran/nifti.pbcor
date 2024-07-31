nifti.pbcor <-
function (nifti1, nifti2, mask_coords, n.parcels = 90, n.parcellations = 50, 
    kmeans.iter.max = 30, kmeans.nstart = 1, kmeans.algorithm = "Hartigan-Wong", 
    kmeans.trace = FALSE, cor.alternative = "two.sided", cor.method = "pearson", 
    cor.exact = NULL, cor.conf.level = 0.94999999999999996, cor.continuity = FALSE, 
    verbose = TRUE) 
{
    if (!inherits(nifti1, "nifti") || nifti1@dim_[1] != 3) {
        stop("'nifti1' must be an object of class 'nifti' with 3D")
    }
    if (!inherits(nifti2, "nifti") || nifti2@dim_[1] != 3) {
        stop("'nifti1' must be an object of class 'nifti' with 3D")
    }
    if (any(dim(nifti2) != dim(nifti1))) {
        stop("'nifti2' must have the same dimensions as 'nifti1'")
    }
    if (!is.matrix(mask_coords) || !is.numeric(mask_coords) || 
        ncol(mask_coords) != 3 || any(mask_coords < 1) || any(apply(mask_coords, 
        2, max) > dim(nifti1))) {
        stop("'mask_coords' must be a matrix of 3D coordinates, as created by 'nifti.pbcor_mask'")
    }
    if (!is.numeric(n.parcels) || n.parcels < 3) {
        stop("'n.parcels' must be a number >= 3 (but many more recommended)")
    }
    if (!is.numeric(n.parcellations) || n.parcellations < 1) {
        stop("'n.parcellations' must be a number >= 1 (but some more recommended)")
    }
    if (verbose) {
        cat("Parcellating the brain and correlating the parcellations")
    }
    r = list()
    for (i_parcellation in 1:n.parcellations) {
        if (verbose) {
            cat(".")
        }
        mask_coords_parcel_ind = kmeans(mask_coords, n.parcels, 
            iter.max = kmeans.iter.max, nstart = kmeans.nstart, 
            algorithm = kmeans.algorithm, trace = kmeans.trace)$cluster
        parcels_x1 = c()
        parcels_x2 = c()
        for (i_parcel in 1:n.parcels) {
            parcel_coords = mask_coords[which(mask_coords_parcel_ind == 
                i_parcel), ]
            parcel_voxels_x1 = c()
            parcel_voxels_x2 = c()
            for (i in 1:nrow(parcel_coords)) {
                parcel_voxels_x1 = c(parcel_voxels_x1, nifti1[parcel_coords[i, 
                  1], parcel_coords[i, 2], parcel_coords[i, 3]])
                parcel_voxels_x2 = c(parcel_voxels_x2, nifti2[parcel_coords[i, 
                  1], parcel_coords[i, 2], parcel_coords[i, 3]])
            }
            parcels_x1 = c(parcels_x1, mean(parcel_voxels_x1))
            parcels_x2 = c(parcels_x2, mean(parcel_voxels_x2))
        }
        r[[i_parcellation]] = cor.test(parcels_x1, parcels_x2, 
            alternative = cor.alternative, method = cor.method, 
            exact = cor.exact, conf.level = cor.conf.level, continuity = cor.continuity)
    }
    if (verbose) {
        cat("\n")
    }
    restimates = sapply(r, function(x) {
        x$estimate
    })
    median_r = r[[which.min(abs(restimates - median(restimates)))]]
    attr(median_r, "parcellations.cor.test") = r
    median_r
}
nifti.pbcor_mask <-
function (nifti, verbose = TRUE) 
{
    if (!inherits(nifti, "nifti") || nifti@dim_[1] != 3) {
        stop("'nifti' must be an object of class 'nifti' with 3D")
    }
    if (verbose) {
        cat("Converting the mask into a coordinates matrix")
    }
    mask_coords = NULL
    for (i in 1:dim(nifti)[1]) {
        if (verbose & i%%3 == 0) {
            cat(".")
        }
        for (j in 1:dim(nifti)[2]) {
            for (k in 1:dim(nifti)[3]) {
                if (nifti[i, j, k] == 1) {
                  mask_coords = rbind(mask_coords, c(i, j, k))
                }
            }
        }
    }
    if (verbose) {
        cat("\n")
    }
    mask_coords
}
