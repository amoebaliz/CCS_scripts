library(mgcv)
library(ncdf4)


# fxns
    replace.fillvalues <- function(arr, ncfile, dimname)
    {
        # change fill value from the ncdf file to NA in R
        fillvalue <- ncatt_get(ncfile, dimname, "_FillValue")$value
        arr[arr == fillvalue] <- NA
        arr
    }


    zmatrix.to.cols <- function (x, y=NULL, zmat=NULL, remove.nas=TRUE) 
    {
        # convert matrix to dataframe. 
        
        # I think there must be a standard more efficient function
        # than this somewhere but I haven't seen it. This is what I
        # have been using for years.

        nmes <- c("x", "y", "z")
        if (is.matrix(x)) 
            x <- data.frame(x)
        if (is.list(x)) {
            if (length(intersect(c("x", "y", "z"), names(x))) != 
                3) {
                nmes <- names(x)[1:3]
                names(x)[1:3] <- c("x", "y", "z")
            }
            zmat <- x$z
            y <- x$y
            x <- x$x
        }
        dat <- expand.grid(x = x, y = y)
        dat$z <- as.vector(zmat)
        if (remove.nas) 
            dat <- dat[!is.na(dat$z), ]
        dat <- dat[order(dat$x, dat$y), ]
        row.names(dat) <- 1:nrow(dat)
        names(dat) <- nmes
        dat
    }


    cols.to.zmatrix <- function (x, y = NULL, z = NULL) 
    {
        # df -> matrix, opposite of zmatrix.to.cols
        nmes <- c("x", "y", "z")
        if (is.matrix(x)) 
            x <- data.frame(x)
        if (is.data.frame(x)) {
            if (length(intersect(c("x", "y", "z"), names(x))) != 
                3) {
                nmes <- names(x)[1:3]
                names(x)[1:3] <- c("x", "y", "z")
            }
            z <- x$z
            y <- x$y
            x <- x$x
        }
        if (any(duplicated(cbind(x, y)))) {
            stop("x-y positions are not unique")
        }
        if ((length(x) != length(y)) | (length(x) != length(z)) | 
            (length(y) != length(z))) {
            stop("x, y and z must be the same length")
        }
        i <- as.numeric(as.factor(x))
        j <- as.numeric(as.factor(y))
        x <- sort(unique(x))
        y <- sort(unique(y))
        mat <- matrix(NA, length(x), length(y))
        for (k in 1:length(z)) {
            mat[i[k], j[k]] <- z[k]
        }
        lst <- list(x, y, mat)
        names(lst) <- nmes
        lst
    }


# read data
    sstfile <- nc_open('sst.nc')
    chlafile <- nc_open('chla.nc')
    ekefile <- nc_open('eke.nc')

    # your dimensions and names probably are different, but this is the
    # idea

    # not sure if you have multiple time slabs in your ncdf file, so I
    # used a variable to make it obvious when you're selecting a time.
    # You could turn this into a function with an index argument if
    # you're running a bunch of times
    dteindex <- 1
    dte <- ncvar_get(sstfile, 'time')[dteindex]

    sst <- ncvar_get(sstfile, 'sst', start=c(1, 1, dteindex),
                     count=c(-1, -1, 1))
    logchla <- log(ncvar_get(chlafile, 'chlorophyll',
                             start=c(1, 1, dteindex),
                             count=c(-1, -1, 1)))
    eke <- ncvar_get(ekefile, 'eke', start=c(1, 1, dteindex),
                     count=c(-1, -1, 1))

    sst <- replace.fillvalues(sst, sstfile, 'sst')
    logchla <- replace.fillvalues(logchla, chlafile, 'chlorophyll')
    eke <- replace.fillvalues(eke, ekefile, 'eke')
    
    # note I assumed dimensions and pixel locations are the same for
    # all three predictor variables
    lons <- ncvar_get(chlafile, 'longitude')
    lats <- ncvar_get(chlafile, 'latitude')
    # check the data look ok
	image(lons, lats, logchla, col=terrain.colors(10))
	image(lons, lats, sst, col=rev(heat.colors(10)))
	image(lons, lats, eke, col=rev(heat.colors(10)))
           
    dat <- zmatrix.to.cols(lons, lats, sst, remove.nas=F)
    names(dat) <- c('lon', 'lat', 'sst')
    dat$logchla <- zmatrix.to.cols(lons, lats, logchla, remove.nas=F)$z
    dat$eke <- zmatrix.to.cols(lons, lats, eke, remove.nas=F)$z
    
    nc_close(sstfile)
    nc_close(chlafile)
    nc_close(ekefile)
 

# predict on Karen's model
    load('cufesmodel.rda')
    # the model was saved as 'cufesmodel'

    # there's a small hitch in using Karen's model that I had
    # forgotten about. We blocked by year to get at the eke/eddy stuff
    # a little better. I suggest using an 'average' year in her model.
    # the median year in her developmental data was 2004, and the
    # average was ~2003.8. So I'd go with 2004.
    dat$year <- 2004
    dat$year <- factor(dat$year, levels <- 1998:2009)
    dat$preds <- plogis(predict(cufesmodel, newdata=dat))


# If you want to go back to ncdf with predictions:
    # output in matrix format
    probs <- cols.to.zmatrix(dat[c('lon', 'lat', 'preds')])$preds
    image(probs)

    # create a new netcdf file
    missval <- -1
    timedim <- ncdim_def('time',
                         units='seconds since 1970-01-01T00:00:00Z',
                         vals=dte,
                         longname='Centered Time', calendar='gregorian')
    latdim <- ncdim_def('latitude', units='degrees_north',
                        vals=ncvar_get(chlafile, 'latitude'),
                        longname='Latitude')
    londim <- ncdim_def('longitude', units='degrees_east',
                        vals=ncvar_get(chlafile, 'longitude'),
                        longname='Longitude')
    preds <- ncvar_def('prediction', units='probability, 0 to 1', 
                        dim=list(timedim, londim, latdim),
                        missval=missval,
                        longname='Predicted probabilities of capture',
                        compression=9, chunksizes=c(1, length(londim$vals),
                                                    length(latdim$vals)))
    predfile <- nc_create('sardine_predictions.nc', preds, force=T)
    ncvar_put(predfile, preds, probs, start=c(dteindex, 1, 1),
              count=c(1, -1, -1), verbose=F)
 
    # attributes aren't strictly needed but if you want to get fancy...
    ncatt_put(predfile, varid=0, 'description',
              paste('Predicted probabilities of capturing',
                    'sardine eggs from Nieto et al model'))
    
    ncatt_put(predfile, 'latitude', '_CoordinateAxisType', 'Lat')
    ncatt_put(predfile, 'latitude', 'coordsys', 'geographic')
    ncatt_put(predfile, 'latitude', 'standard_name', 'latitude')
    ncatt_put(predfile, 'latitude', 'axis', 'Y')
    ncatt_put(predfile, 'latitude', 'spacing', 'even')

    ncatt_put(predfile, 'longitude', '_CoordinateAxisType', 'Lon')
    ncatt_put(predfile, 'longitude', 'coordsys', 'geographic')
    ncatt_put(predfile, 'longitude', 'standard_name', 'longitude')
    ncatt_put(predfile, 'longitude', 'axis', 'X')
    ncatt_put(predfile, 'longitude', 'spacing', 'even')

    ncatt_put(predfile, 'time', '_CoordinateAxisType', 'Time')
    ncatt_put(predfile, 'time', 'standard_name', 'time')
    ncatt_put(predfile, 'time', 'axis', 'T')

    ncatt_put(predfile, 'prediction', 'colorBarMaximum', 1)
    ncatt_put(predfile, 'prediction', 'colorBarMinimum', 0)

    nc_close(predfile)

