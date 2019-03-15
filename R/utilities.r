# suppressPackageStartupMessages(require(zeallot))
# suppressPackageStartupMessages(require(lumi))
# suppressPackageStartupMessages(require(methylumi))
# suppressPackageStartupMessages(require(IlluminaHumanMethylation450kprobe))
# suppressPackageStartupMessages(require("FDb.InfiniumMethylation.hg19"))
# suppressPackageStartupMessages(require(IlluminaHumanMethylationEPICmanifest))
# suppressPackageStartupMessages(require("IlluminaHumanMethylationEPICanno.ilm10b2.hg19"))
# suppressPackageStartupMessages(require("sva"))
# suppressPackageStartupMessages(require("wateRmelon"))
# suppressPackageStartupMessages(require("parallel"))
# suppressPackageStartupMessages(require("ggplot2"))
# suppressPackageStartupMessages(require("logging"))
# suppressPackageStartupMessages(require(digest))
# suppressPackageStartupMessages(require(RColorBrewer))
# suppressPackageStartupMessages(require(heatmap3))

formatSexNames = function(sex_vector){
    ## in F/M/U format
    new_vector = NULL
    for (s in sex_vector){
        if (substr(s,1,1) %in% c('F', 'f')){
            new_vector = c(new_vector, 'F')
        } else if (substr(s,1,1) %in% c('M', 'm')){
            new_vector = c(new_vector, 'M')
        } else {
            new_vector = c(new_vector, 'U')
        }
    }
    return (new_vector)
}
##

getBarcodes = function(filenames, suffix=c('_Grn.idat', '_Red.idat')){
    ## extract the unique alnum names; if the pos is not continuous, include
    ## the part(s) between
    # return(unique(sapply(filenames,
    return(unique(vapply(filenames,
                         function(x) sub(paste(suffix, collapse='|'),
                                         '',
                                         basename(x)), "")))
}

##
getUniqueNames = function(name_vectors){
    ## extract the unique alnum names; if the pos is not continuous, include the part(s) between
    char_mat = do.call(rbind,
                       lapply(name_vectors,
                              function(x) {
                                  unlist(strsplit(x,'[^[:alnum:]]+'))
                              }))
    sep_mat = unique(do.call(rbind,
                             lapply(name_vectors,
                                    function(x) {
                                        unlist(strsplit(x,'[[:alnum:]]+'))
                                    })))
    if (nrow(sep_mat)>1) {
        sep_vector = rep('_', ncol(char_mat))
    } else {
        sep_vector = sep_mat[1,]
    }
    slen = apply(char_mat, 2, function(x) {length(unique(x))})
    column_list = list()
    pasted_vector = NULL
    for (i in order(slen, decreasing=TRUE)){
        if (!length(column_list)) {
            column_list[[as.character(i)]] = char_mat[,i]
            pasted_vector = char_mat[,i]
        } else {
            pasted_vector = paste(pasted_vector, char_mat[,i], sep='::')
            column_list[[as.character(i)]] = char_mat[,i]
        }
        if (length(unique(pasted_vector)) == length(name_vectors)){
            break
        }
    }

    pasted_vector = NULL
    # sorted_cols = sort(sapply(names(column_list),
    sorted_cols = sort(vapply(names(column_list),
                              function(x) as.integer(x), ""))
    for (i in sorted_cols) {
        if(is.null(pasted_vector)) {
            pasted_vector = column_list[[as.character(i)]]
        } else {
            pasted_vector = paste(pasted_vector,
                                  column_list[[as.character(i)]],
                                  sep=sep_vector[i])
        }
    }
    return(pasted_vector)
}

##
#' @importFrom Biobase fData fData<- featureNames
#' @importFrom AnnotationDbi mget keys
#' @importFrom GenomicFeatures features
#' @importFrom GenomeInfoDb seqnames
#' @importFrom BiocGenerics start
addAnnotation = function (methyLumiM, lib = "FDb.InfiniumMethylation.hg19",
                          annotationColumn = c("COLOR_CHANNEL",
                                               "CHROMOSOME",
                                               "POSITION",
                                               "CN_TYPE")) {
    ### modify lumi::addAnnotationInfo to add "CN_TYPE"
    if (is(methyLumiM, "MethyLumiM") | is(methyLumiM, "MethyLumiSet")) {
        ff <- fData(methyLumiM)
        probeList <- featureNames(methyLumiM)
    }else if (is.character(methyLumiM)) {
        probeList <- methyLumiM
        ff <- data.frame()
    } else {
        stop("methyLumiM should be either a MethyLumiM object of a character vectors of probe names!")
    }
    # if (!is.null(lib) && require(lib, character.only = TRUE)) {
    if (!is.null(lib)) {
        ff <- NULL
        if (exists(lib)) {
            lib <- get(lib)
            if (is(lib, "FeatureDb")) {
                allAnnotation <- features(lib)
            }
            ff <- data.frame(ProbeID = probeList,
                             CHROMOSOME = NA,
                             POSITION = NA, COLOR_CHANNEL = NA,
                             CN_TYPE = NA)
            rownames(ff) <- probeList
            if (any(!(probeList %in% names(allAnnotation)))) {
                missingProbe <- probeList[!(probeList %in% names(allAnnotation))]
                probeList <- probeList[probeList %in% names(allAnnotation)]
                warnings(paste(paste(missingProbe, collapse = ","),
                               "probes do not exist in the annotation library!"))
            }
            allAnnotation <- allAnnotation[probeList]
            ff[probeList, "CHROMOSOME"] <- as.character(seqnames(allAnnotation))
            ff[probeList, "POSITION"] <- as.numeric(start(allAnnotation))
            ff[probeList, "COLOR_CHANNEL"] <- as.character(allAnnotation$channel450)
            ff[probeList, "CN_TYPE"] <- as.character(allAnnotation$probeType)
        }
        if (is.null(ff)) {
            colorInfo <- chr <- loc <- rep(NA, length(probeList))
            names(colorInfo) <- names(chr) <- names(loc) <- probeList
            obj <- get(paste(sub("\\.db$", "", lib), "COLORCHANNEL",
                             sep = ""))
            pp <- probeList[probeList %in% keys(obj)]
            # colorInfo[pp] <- sapply(AnnotationDbi::mget(pp, obj),
            colorInfo[pp] <- vapply(mget(pp, obj), function(x) x[1], "")
            ff$COLOR_CHANNEL <- colorInfo
            chrPattern <- paste(sub("\\.db$", "", lib), "CHR", sep = "")
            chrObjs <- ls(paste("package:", lib, sep = ""),
                          pattern = paste(chrPattern,
                                          "[0-9]+$", sep = ""))
            if (length(chrObjs) > 1) {
                chrVersion <- sub(".*[^0-9]([0-9]+)$", "\\1", chrObjs)
                obj <- get(chrObjs[which.max(as.numeric(chrVersion))])
            } else {
                obj <- get(chrPattern)
            }
            pp <- probeList[probeList %in% keys(obj)]
            # chr[pp] <- sapply(AnnotationDbi::mget(pp, obj), function(x) x[1])
            chr[pp] <- vapply(mget(pp, obj), function(x) x[1], "")
            ff$CHROMOSOME <- as.character(chr)
            obj <- get(paste(sub("\\.db$", "", lib), "CPGCOORDINATE",
                             sep = ""))
            pp <- probeList[probeList %in% keys(obj)]
            # loc[pp] <- sapply(AnnotationDbi::mget(pp, obj), function(x) x[1])
            loc[pp] <- vapply(mget(pp, obj), function(x) x[1], "")
            ff$POSITION <- as.numeric(as.character(loc))
        }
    } else {
        if (all(c("CHR", "MAPINFO") %in% names(ff))) {
            ff$CHROMOSOME <- ff$CHR
            ff$POSITION <- as.numeric(ff$MAPINFO)
        }
        if (!all(c("CHROMOSOME", "POSITION", "COLOR_CHANNEL", "CN_TYPE") %in%
                 names(ff))) {
            stop("Probe annotation information is not available. Please provide annotation library!")
        }
    }
    if (!is.data.frame(ff))
        ff <- as.data.frame(ff, stringsAsFactors = FALSE)
    rownames(ff) <- ff$ProbeID
    if (is(methyLumiM, "MethyLumiM")) {
        fData(methyLumiM) <- ff
        return(methyLumiM)
    }else {
        return(ff)
    }
}


##
#' @importFrom stats dendrapply as.dendrogram is.leaf
#' @importFrom graphics plot
colorCluster = function(hc, labs=NULL, lab.cols=NULL, tip.shapes=NA,
                        show.labs=TRUE, main.txt='') {

    #mult.tips = FALSE
    if (is.null(tip.shapes)){ tip.shapes=NA }
    if (length(labs)>0){
        if (length(labs) != length(lab.cols)) {
            stop('--lab(el)s and lab(el)col(or)s must match in numbers')
        }
    }
    if (is.null(labs)){
        if (!identical(names(lab.cols), hc$labels)) {
            stop('--if no labs provided, vector lab.cols must have the sample names as its names')
        }
        if (length(tip.shapes)>1){
            #mult.tips = TRUE
            if (!identical(names(tip.shapes), hc$labels)) {
                stop('--if no labs provided, vector tip.shapes must have the sample names as its names')
            }
        }
    }

    dList <- dendrapply(as.dendrogram(hc), function(n) {
        if(is.leaf(n)){
            # print(attributes(n)$nodePar)
            if (is.null(labs)){
                # labelCol <- paste("#",substring(digest(attr(n,"label")),1,6), sep="")
                labelCol <- lab.cols[attr(n,"label")] #
                if (length(tip.shapes)>1){
                    tip <- tip.shapes[attr(n,"label")] #
                } else {
                    tip <- tip.shapes
                }
            } else {
                labelCol <- lab.cols[which(labs == attr(n,"label"))]
                if (length(tip.shapes)>1) {
                    tip <- tip.shapes[which(labs == attr(n,"label"))] #
                } else {
                    tip <- tip.shapes
                }
            }
            attr(n, "edgePar") <- list(col = labelCol)
            if (show.labs){ # how to color the tip shapes???
                attr(n, "nodePar") <- list(pch = tip,
                                           col = labelCol,
                                           lab.col = labelCol,
                                           lab.cex = 0.50)
            } else {
                attr(n, "nodePar") <- list(pch = tip,
                                           col = labelCol,
                                           lab.col = 'white',
                                           lab.cex = 0.0001)
            }
        }
        n;
    })
    plot(dList, xlab='', ylab='distance', main=main.txt)
}

##
#plotClustChrX = plotQC(mset.lumi[fData(mset.lumi)$CHROMOSOME=='chrX',], groupCol='Sex', plots=c('mds', 'hclust'),outdir, outname='chrX')

##
sampleByGroup = function(sampleVector, groupVector, maxNum){
    if(length(sampleVector) > maxNum && length(unique(groupVector)) > 1) {
        counts = sort(table(groupVector))
        #counts = as.numeric(tmp)
        #names(counts) = names(tmp)
        allsubrows = NULL
        for ( i in seq_len(length(counts))){
            grp = names(counts)[i]
            val = counts[grp]
            gnum = ceiling((maxNum-length(allsubrows)) / (length(unique(groupVector))-i+1))
            num = ifelse(val > gnum,
                         gnum,
                         val)
            #subsamples = which(groupVector==grp)[1:num]
            subrows = sample(which(groupVector==grp),
                             num,
                             replace=FALSE)
            allsubrows = c(allsubrows, subrows[!is.na(subrows)])
        }
        return(sampleVector[allsubrows])
    } else if (length(sampleVector) > maxNum) {
        return (sample(sampleVector,
                       maxNum,
                       replace=FALSE))
    } else {
        return (sampleVector)
    }
}

## use this function to plot: color bias plot, density plot, and lumi-outlier
#' @importFrom lumi boxplotColorBias
#' @importFrom grDevices pdf dev.off
#' @importFrom stats quantile cmdscale hclust cor dist median as.dist cutree
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics plot legend abline legend
#' @importFrom BiocGenerics density
plotQC = function(mset.lumi, outdir='.', step=NULL, outname=NULL, scaled=TRUE,
                  plots=c('color',
                          'density',
                          'outlier'),
                  maxnum.samples=NULL, sampling.variable=NULL,
                  cluster.variable=NULL,
                  distance.method=c('correlation',
                                    "euclidean",
                                    "maximum",
                                    "manhattan",
                                    "canberra",
                                    "binary",
                                    "minkowski"),
                  cluster.method=c("average",
                                   "ward.D2",
                                   "single",
                                   "complete",
                                   "mcquitty",
                                   "median",
                                   "centroid"),
                  outlier.threshold=2, outlier.plot=TRUE) {
    stopifnot(is(mset.lumi, "ExpressionSet"))
    list2return = list()

    distance.method = match.arg(distance.method)
    cluster.method = match.arg(cluster.method)

    ## outlier detection will use all samples/columns
    if ('outlier' %in% plots){
        mdat = exprs(mset.lumi)
    }
    ## set colors by group
    grp_colors=NULL
    if (!is.null(cluster.variable)){
        if (length(cluster.variable)>1){
            cluster.variable = cluster.variable[1]
            warning(sprintf('***** only the first one of [%s] is selected as the group name',
                            paste(cluster.variable, collapse=',')))
        }
        grp_levels = as.character(sort(unique(pData(mset.lumi)[[cluster.variable]])))
        grp_vector = as.character(pData(mset.lumi)[[cluster.variable]])
        grp_colors = brewer.pal(n = 8, name = "Dark2")
        if (length(grp_levels)>8){
            grp_colors = rep(grp_colors, ceiling(length(grp_levels)/8))
        }
        grp_colors = grp_colors[seq_len(length(grp_levels))]
        names(grp_colors) = grp_levels
        col_vector = grp_colors[grp_vector]; names(col_vector) = sampleNames(mset.lumi)
    } else {
        col_vector = rep('black', length(sampleNames(mset.lumi))); names(col_vector) = sampleNames(mset.lumi)
        grp_colors = NULL
    }
    ## subset columns if needed
    list2return[['subset.samples']] = NULL
    if (!is.null(maxnum.samples)){
        if(maxnum.samples >= length(sampleNames(mset.lumi))){
            warning(sprintf('***** only [%s] samples in data', length(sampleNames(mset.lumi))))
        } else {
            if (is.null(sampling.variable) && is.null(cluster.variable)){
                subcols = sample(seq_len(length(sampleNames(mset.lumi))),
                                 maxnum.samples,
                                 replace=FALSE)
            } else {
                if (!is.null(cluster.variable)){
                    sampling.variable = cluster.variable
                }
                subcols = sampleByGroup(seq_len(length(sampleNames(mset.lumi))),
                                        pData(mset.lumi)[[sampling.variable]],
                                        maxnum.samples)
            }
            mset.lumi = mset.lumi[,subcols]
            list2return[['subset.samples']]= sampleNames(mset.lumi)
        }
    }

    if ('color' %in% plots){ ### need to change to using ggplot2!!!!
        print('-----color bias boxplot')
        for (channel in c('unmethy', 'methy', 'both', 'sum')){
            pdf(sprintf ('%s/QC%s%s_%s_colorBias.pdf',
                         outdir,
                         step,
                         ifelse(is.null(outname),
                                '',
                                paste('_', outname, sep='')),
                         channel),
                width=5,
                height=5)
            ## if want the R/G values, use boxplotColorBias_re to return, o.w. use original function boxplotColorBias
            bias.RG = lumi::boxplotColorBias(mset.lumi,
                                       channel=channel,
                                       main=sprintf('Color bias for %s channels',
                                                    channel))
            dev.off()
        }
    }

    if ('density' %in% plots){
        print ('-----density plot')
        pdf(sprintf ('%s/QC%s%s_density.pdf',
                     outdir,
                     step,
                     ifelse(is.null(outname),
                            '',
                            paste('_', outname, sep=''))))

        dens = BiocGenerics::density(mset.lumi,
                              na.rm = TRUE,
                              xlab="M-Value",
                              main='Density of M-values',
                              legend=ifelse(length(sampleNames(mset.lumi))>=30,
                                            c(-1000,1000), 'right'))
        dev.off()
    }
    ####
    if ('outlier' %in% plots){
        #### modified from the function detectOutlier in package lumi
        #Th = outlier.threshold #, ifPlot = FALSE
        # create a 'center' pseudo-sample using rowmeans
        # calc dist b/w 'center' and each sample
        # get the median of these dist or mad.d
        # exclude 10% samples of the largest dist
        # update the 'center' pseudo-sample
        # update mad.d (the median of those center-vs-sample dist)
        # set a threshold of pre-set "outlier.threshold * mad.d"
        # any sample of which the dist to center is larger than this threshold is deemed as an outlier
        # I'd use clustering to further the selection
        print('-----outlier detection based on clustering')
        center <- rowMeans(mdat, na.rm = TRUE)
        mdat_add_center <- cbind(center, mdat)
        colnames(mdat_add_center) <- c("Center", colnames(mdat))
        if (grepl(sprintf('^%s',distance.method), 'correlation',
                  ignore.case = TRUE)){
            d <- 1 - cor(mdat_add_center,
                         use="pairwise.complete.obs")
            #mad.d <- median(1 - d[2:nrow(d), 1])
        } else {
            if (scaled) {
                mdat_add_center <- scale(mdat_add_center)
            }
            d <- as.matrix(dist(t(mdat_add_center), method = distance.method))
            #mad.d <- median(d[2:nrow(d), 1])
        }
        mad.d <- median(d[2:nrow(d), 1])
        excludeInd <- which(rank(d[2:nrow(d), 1]) > ncol(mdat) * 0.99)
        # d is transposed
        center <- rowMeans(mdat[, -excludeInd],
                           na.rm = TRUE) # update center
        mdat_add_center <- cbind(center, mdat) # update mdat_add_center
        colnames(mdat_add_center) <- c("Center", colnames(mdat))
        if (grepl(sprintf('^%s',distance.method),
                  'correlation',
                  ignore.case = TRUE)) {
            d <- 1- cor(mdat_add_center, use="pairwise.complete.obs")
        } else {
            if (scaled) {
                mdat_add_center <- scale(mdat_add_center)
            }
            d <- as.matrix(dist(t(mdat_add_center), method = distance.method))
        }
        mad.d <- median(d[2:nrow(d), 1])  # update mad.d
        outlier_labels <- (d[2:nrow(d), 1] >= outlier.threshold * mad.d)
        dist_outlier_names = names(outlier_labels)[outlier_labels==TRUE]

        hc = hclust(as.dist(d),
                    method=cluster.method) ##
        list2return$dist.hclust = hc
        tip_shapes = rep(1, ncol(mdat_add_center))
        names(tip_shapes) = colnames(mdat_add_center)
        if (sum(outlier_labels) > 0){
            hc_labels = cutree(hc, 2)
            stopifnot(identical(names(outlier_labels), names(hc_labels[2:length(hc_labels)])))
            hc_labels[hc_labels==hc_labels[1]] = FALSE
            hc_labels[hc_labels!=hc_labels[1]] = TRUE
            #hc_labels = as.logical(hc_labels)

            tmpdf = data.frame(distThreshold=outlier_labels,
                               distCluster=as.logical(hc_labels[2:length(hc_labels)]),
                               outlier=FALSE)
            # the 1st is the 'Center', which must be a non-outlier
            clust_outlier_names = names(hc_labels)[hc_labels!=hc_labels[1]]
            if (identical(setdiff(dist_outlier_names,clust_outlier_names),  character(0))){
                dist_outlier_names = clust_outlier_names
            }

            tmpdf[dist_outlier_names, 'outlier'] =TRUE
            #attr(outlier, "sampleDistance") <- d
            #attr(outlier, "threshold") <- outlier.threshold #??

            tip_shapes[dist_outlier_names] = 16

            list2return$dist.outliers = tmpdf
        } else {
            list2return$dist.outliers = NULL
        }

        main <- paste("Outlier detection based on sample distance to \"Center\"")
        pdf(sprintf ('%s/QC%s%s_outliers.pdf',
                     outdir,
                     step,
                     ifelse(is.null(outname),
                            '',
                            paste('_', outname, sep=''))),
            width=10,
            height=8)
        colorCluster(hc,
                     main.txt=main,
                     lab.cols=c('Center'='gray', col_vector),
                     tip.shapes=tip_shapes)
        abline(h = outlier.threshold * mad.d, col = 'red', lty = 2)
        if (!is.null(grp_colors)){
            legend('top', names(grp_colors),
                   text.col=grp_colors,
                   bty='o',
                   bg='white',
                   cex=1.0,
                   ncol=length(grp_colors))
        }
        dev.off()
    }

    ####
    return(list2return)
}

### heatmap3 of mset.lumi obj
#' @importFrom stats as.dist cor dist
distFun = function(x) {
    if (grepl('^cor', x, ignore.case = TRUE)){
        return (function(m)(as.dist(1 - cor(t(m), use = "pa"))))
    } else if (sum(grepl(x, c("euclidean",
                              "maximum",
                              "manhattan",
                              "canberra",
                              "binary",
                              "minkowski")))==1){
        return (function(m)(dist(m, method = x, diag = FALSE, upper = FALSE, p = 2)))
    } else {
        stop('----- function name must be a unique frontal substring of "cor", "euclidean", "maximum", "manhattan", "canberra", "binary", and "minkowski"')
    }
}

## issues: remove label txt; remove legend plot, add legend for colors
# groupByChrom = plotHeatmap(mset.lumi[fData(mset.lumi)$CHROMOSOME==chrom,],
#                            outname=chrom, plots=c('heatmap', 'mds', 'hclust'),
#                            maxnum.samples=NULL,
#                            sampling.variable=NULL,
#                            sidebar.variables=sexColName,
#                            cluster.variable=sexColName,
#                            distance.method="euclidean",
#                            feature.selection.method = 'mad',
#                            feature.selection.quantile = NULL,
#                            scaled='row',
#                            cluster.method="ward.D2",
#                            balanceColor=TRUE,
#                            cexCol=1,
#                            cexRow=1,
#                            col.margin=5,
#                            row.margin=5,
#                            labRow=NA,
#                            labCol=NA)

#' @importFrom Biobase pData sampleNames exprs
#' @importFrom AnnotationDbi sample
#' @importFrom heatmap3 heatmap3
#' @importFrom stats quantile cmdscale hclust
#' @importFrom graphics plot legend
#' @importFrom grDevices dev.off
plotHeatmap_lumi = function(mset.lumi, outname=NULL,
                            plots=c('heatmap', 'mds', 'hclust'),
                            maxnum.samples=NULL, sampling.variable=NULL,
                            sidebar.variables=NULL, cluster.variable=NULL,
                            distance.method=c("euclidean",
                                              "maximum",
                                              "manhattan",
                                              "canberra",
                                              "binary",
                                              "minkowski",
                                              'correlation'),
                            feature.selection.method = c('mad', 'sd'),
                            feature.selection.quantile = NULL,
                            scaled=c('row', 'column', 'none'),
                            cluster.method=c("ward.D2",
                                             "single",
                                             "complete",
                                             "average",
                                             "mcquitty",
                                             "median",
                                             "centroid"),
                            balanceColor=TRUE,
                            cexCol=1, cexRow=1, col.margin=5, row.margin=5,
                            labRow=c(NA, NULL), labCol=c(NA, NULL),
                            highlights=NULL){
    # require(RColorBrewer)
    stopifnot(is(mset.lumi, "ExpressionSet"))
    distance.method = match.arg(distance.method)
    cluster.method = match.arg(cluster.method)
    feature.selection.method = match.arg(feature.selection.method)
    scaled = match.arg(scaled)
    #cluster.info = match.arg(clust.info)

    if (is.null(sampling.variable) &&
        is.null(cluster.variable) &&
        ('mds' %in% plots || 'hclust' %in% plots)){
        stop(sprintf('***** either cluster.variable or sampling.variable must be specified to do clustering'))
    }

    stopifnot(feature.selection.quantile <= 1 || is.null(feature.selection.quantile))
    if (!is.null(feature.selection.quantile) && feature.selection.quantile<0.5){
        warning(sprintf('----provided quantile is %s that is less than 0.5; will select top features more variable than %s of the all',
                        feature.selection.quantile, feature.selection.quantile))
    }
    list2return = list()

    ## check if the colnames in para exist(s)
    check_colnames = NULL
    for (variable in unique(c(sidebar.variables,sampling.variable,cluster.variable))){
        if (!variable %in% colnames(pData(mset.lumi))){
            check_colnames = c(check_colnames, variable)
        }
    }
    if (!is.null(check_colnames)){
        stop(sprintf('***** [%s] is/are not in pheno data %s',
                     paste(check_colnames, collapse=','),
                     deparse(substitute(mset.lumi))))
    }

    ## subset columns if needed
    if (!is.null(maxnum.samples) && maxnum.samples < length(sampleNames(mset.lumi))){
        if (is.null(sampling.variable) && is.null(cluster.variable)){
            subcols = sample(seq_len(length(sampleNames(mset.lumi))),
                             maxnum.samples, replace=FALSE)
        } else {
            if (!is.null(cluster.variable)){
                sampling.variable = cluster.variable
            }
            subcols = sampleByGroup(seq_len(length(sampleNames(mset.lumi))),
                                    pData(mset.lumi)[[sampling.variable]],
                                    maxnum.samples)
        }
        mset.lumi = mset.lumi[,subcols]
    }

    mdat = exprs(mset.lumi)
    #feature.selection.quantile = 0.999
    ## subset rows if needed
    if (!is.null(feature.selection.quantile)){
        tmp <- apply(mdat, 1, get(feature.selection.method), na.rm=TRUE)
        subrows = which(tmp>=quantile(tmp, feature.selection.quantile))
        if (length(subrows)<2){
            stop(sprintf('----- too few rows remained after filtering for top %s most variable rows',
                         feature.selection.quantile))
        }
        mset.lumi = mset.lumi[subrows,]
        mdat = mdat[subrows,]
    }

    ## group colors for hclust and/or mds plot
    if ('mds' %in% plots || 'hclust' %in% plots){
        mdist = distFun(distance.method)(t(mdat))
        grp_colors=NULL
        if (!is.null(cluster.variable)){
            if (length(cluster.variable) > 1){
                cluster.variable = cluster.variable[1]
                warning(sprintf('***** only the first one of [%s] is selected as the group name',
                                paste(cluster.variable, collapse=',')))
            }
            grp_levels = as.character(sort(unique(pData(mset.lumi)[[cluster.variable]])))
            grp_vector = as.character(pData(mset.lumi)[[cluster.variable]])
            grp_colors = brewer.pal(n = 8, name = "Dark2")
            if (length(grp_levels)>8){
                grp_colors = rep(GroupColors, ceiling(length(grp_levels)/8))
            }
            grp_colors = grp_colors[seq_len(length(grp_levels))]
            names(grp_colors) = grp_levels
            #col_vector = grp_colors[grp_vector]
        }

        if(!is.null(grp_colors)){
            tmpcols = grp_colors[grp_vector]
        } else {
            tmpcols = 'black'
        }

        ## highlight samples e.g. outliers
        shapes = rep(1, ncol(mdat))
        names(shapes) = colnames(mdat)
        if (!is.null(highlights)){
            if (all(highlights %in% colnames(mdat))){
                shapes[highlights] = 16
            } else {
                warning(sprintf('----- samples/assays in highlights are not a subset of sample names of mset.lumi data',
                                paste(highlights, sep=',')))
            }
        }

        if ('mds' %in% plots){
            pdf(sprintf ('%s/QC%s_MDS.pdf',
                         outdir,
                         ifelse(is.null(outname),
                                '',
                                paste('_', outname, sep=''))),
                width=5,
                height=5)
            #print (cmds)
            cmds = cmdscale(mdist, k=2)
            list2return$cmds = cmds
            plot(cmds[,1],
                 cmds[,2],
                 xlab="Dimension 1",
                 ylab="Dimension 2",
                 col=tmpcols,
                 pch=shapes,
                 lwd=1.5,
                 main='MDS of mset.lumi')
            if(!is.null(grp_colors)){
                legend('top',
                       grp_levels,
                       text.col=grp_colors,
                       bty='o',
                       cex=1.0,
                       ncol=length(grp_levels)) #horiz=T)
            }
            dev.off()
        }
        #colorCluster = function(hc, labs=NULL, lab.cols=NULL, tip.shapes=NULL, show.labs=TRUE, main.txt=''){
        if ('hclust' %in% plots){
            pdf(sprintf ('%s/QC%s_hclust.pdf',
                         outdir, ifelse(is.null(outname),
                                        '',
                                        paste('_', outname, sep=''))),
                width=5,
                height=5)
            hclust.obj = hclust(mdist, method =cluster.method)
            list2return$hclust = hclust.obj #$merge
            colorCluster(hclust.obj, labs=sampleNames(mset.lumi),
                         lab.cols=tmpcols,
                         tip.shapes=shapes,
                         show.labs=TRUE,
                         main.txt='clustering of mset.lumi')
            if(!is.null(grp_colors)){
                legend('top', grp_levels,
                       text.col=grp_colors,
                       bty='o',
                       bg='white',
                       cex=1.0,
                       ncol=length(grp_levels)) #horiz=T)
            }
            dev.off()
        }
    }

    ### for heatmap3
    if ('heatmap' %in% plots){
        if (!is.null(sidebar.variables)){
            # set colors
            ColSideColors = NULL
            variableColors = list()
            #GroupLevelColorList = NULL # list
            for ( variable in sidebar.variables){
                grp_levels = sort(unique(as.character(pData(mset.lumi)[[variable]])))
                grp_vector = as.character(pData(mset.lumi)[[variable]])
                ## colors
                grp_colors = brewer.pal(n = 8, name = "Dark2")
                if (length(grp_levels)>8){
                    grp_colors = rep(grp_colors, ceiling(length(grp_levels)/8))
                }
                grp_colors = grp_colors[seq_len(length(grp_levels))]
                names(grp_colors) = grp_levels
                variableColors[[variable]] = grp_colors
                #for ( mycol in brewColors){
                #    GroupLevelColorList[[variable]][[brewColors[mycol]]] = mycol
                #}
                col_vector = grp_colors[grp_vector] #rep('dark grey', nrow(pData(mset.lumi)))
                #for (g in grp_vector){
                #    col_vector[which(pData(mset.lumi)[[variable]]==g)] = grp_colors[g]
                #}
                ColSideColors = cbind(ColSideColors, col_vector)
                colnames(ColSideColors)[ncol(ColSideColors)] = variable
                #rownames(ColSideColors)=NULL
            }
        }

        pdf(sprintf ('%s/QC%s_heatmap.pdf',
                     outdir,
                     ifelse(is.null(outname),
                            '',
                            paste('_', outname, sep=''))),
            width=5,
            height=5)
        if (is.null(ColSideColors)){
            heatmap3(mdat,
                     showRowDendro=FALSE,
                     scale=scaled,
                     distfun=distFun(distance.method),
                     method=cluster.method,
                     balanceColor=balanceColor,
                     cexCol=cexCol,
                     cexRow=cexRow,
                     margins=c(col.margin,row.margin),
                     labRow=labRow,
                     labCol=labCol)
        } else {
            grpnames=NULL; colornames=NULL
            for (variable in names(variableColors)){
                grpnames=c(grpnames, sprintf('(%s)', variable)); colornames=c(colornames, 'black')
                grpnames=c(grpnames, names(variableColors[[variable]])); colornames=c(colornames,variableColors[[variable]])
            }
            heatmap3(mdat,
                     ColSideColors=ColSideColors,
                     showRowDendro=FALSE,
                     scale=scaled,
                     distfun=distFun(distance.method),
                     method=cluster.method,
                     balanceColor=balanceColor,
                     cexCol=cexCol,
                     cexRow=cexRow,
                     margins=c(col.margin,row.margin),
                     labRow=labRow,
                     labCol=labCol)
            legend(0, 1, xpd=TRUE,
                   grpnames,
                   text.col=colornames,
                   bty='n',
                   cex=1.0,
                   ncol=1,
                   adj=c(1,1))
        }
        dev.off()

        #if (cluster.info && !'hclust' %in% plots){ ## current function won't keep the info
        #    list2return$hclust = hclust(distFun(distance.method)(t(mdat)), method =cluster.method)
        #}
    }
    return (list2return)
}


##
reLabelVector = function(list1, list2){
    ### 2-level vector only; change the labels in list2 same as those in list1 so that two better match
    ele1 = sort(unique(list1)); ele2 = sort(unique(list2))
    idlist = list()
    for ( i in seq_len(2)){
        tmp = ifelse(list2==ele2[i], ele1[1], ele1[2])
        idlist[i] = sum(ifelse(tmp==list1, 1, 0))
    }

    if (idlist[[1]]>idlist[[2]]){
        return (ifelse(list2==ele2[1], ele1[1], ele1[2]))
    } else if (idlist[[1]]<idlist[[2]]) {
        return (ifelse(list2==ele2[2], ele1[1], ele1[2]))
    } else {
        return (NULL)
    }
}

##--------
#' @importFrom stats pf
lmPR <- function (lmobj) {
    if (class(lmobj) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(lmobj)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=FALSE)
    attributes(p) <- NULL
    return(c(p, summary(lmobj)$adj.r.squared))
}

##
colMedianCenter <- function(mat){
    return (sweep(mat, 2, apply(mat, 2, median), "-"))
}

## the converted data from this function is different from using exprs() or betas()!!!!
convertBM <- function(mat){
    if(min(mat[,1])>=0){
        ## B-values to M-values
        tmp = which(mat==1)
        if (length(tmp)>0){
            mat[tmp] = 0.99999
        }
        return (log2(mat/(1-mat)))
    } else {
        return (2^mat/(2^mat+1))
    }
}


###new
#' @importFrom Biobase pData sampleNames exprs
#' @importFrom AnnotationDbi sample
#' @importFrom heatmap3 heatmap3
#' @importFrom stats quantile cmdscale hclust
#' @importFrom RColorBrewer brewer.pal
#' @importFrom graphics plot legend
#' @importFrom grDevices dev.off
plotHeatmap = function(mData, pData=NULL, outname=NULL, outdir=NULL,
                       plots=c('heatmap', 'mds', 'hclust'),
                       maxnum.samples=NULL,
                       sampling.variable=NULL, sidebar.variables=NULL,
                       cluster.variable=NULL,
                       distance.method=c("euclidean",
                                         "maximum",
                                         "manhattan",
                                         "canberra",
                                         "binary",
                                         "minkowski",
                                         'correlation'),
                       feature.selection.method = c('mad', 'sd'),
                       feature.selection.quantile = NULL,
                       scaled=c('row', 'column', 'none'),
                       cluster.method=c("ward.D2",
                                        "single",
                                        "complete",
                                        "average",
                                        "mcquitty",
                                        "median",
                                        "centroid"),
                       balanceColor=TRUE, cexCol=1, cexRow=1, col.margin=5,
                       row.margin=5, labRow=c(NA, NULL), labCol=c(NA, NULL),
                       highlights=NULL){
    # require(heatmap3); require(RColorBrewer)
    stopifnot(is(mData, "ExpressionSet")|| !(is.null(pData)))
    distance.method = match.arg(distance.method)
    cluster.method = match.arg(cluster.method)
    feature.selection.method = match.arg(feature.selection.method)
    scaled = match.arg(scaled)
    #cluster.info = match.arg(clust.info)

    if (is.null(pData)){ pData = pData(mData)}
    if (is(mData, "ExpressionSet")) {
        sampleNames = sampleNames(mData)
    } else {
        sampleNames = colnames(mData)
    }

    if (is.null(sampling.variable) && is.null(cluster.variable) && ('mds' %in% plots || 'hclust' %in% plots)){
        stop(sprintf('***** either cluster.variable or sampling.variable must be specified to do clustering'))
    }

    stopifnot(feature.selection.quantile<=1 || is.null(feature.selection.quantile))
    if (!is.null(feature.selection.quantile) && feature.selection.quantile<0.5){
        warning(sprintf('----provided quantile is %s that is less than 0.5; will select top features more variable than %s of the all',
                        feature.selection.quantile,
                        feature.selection.quantile))
    }
    list2return = list()

    ## check if the colnames in para exist(s)
    check_colnames = NULL
    for (variable in unique(c(sidebar.variables,
                              sampling.variable,
                              cluster.variable))){
        if (!variable %in% colnames(pData)){
            check_colnames = c(check_colnames, variable)
        }
    }
    if (!is.null(check_colnames)){
        stop(sprintf('***** [%s] is/are not in pheno data %s',
                     paste(check_colnames, collapse=','),
                     ifelse(is(mData, "ExpressionSet"),
                            deparse(substitute(mData)),
                            deparse(substitute(pData)) ) ))
    }

    ## subset columns if needed
    if (!is.null(maxnum.samples) && maxnum.samples < length(sampleNames)){
        if (is.null(sampling.variable) && is.null(cluster.variable)){
            subcols = sample(seq_len(length(sampleNames)),
                             maxnum.samples, replace=FALSE)
        } else {
            if (!is.null(cluster.variable)){
                sampling.variable = cluster.variable
            }
            subcols = sampleByGroup(seq_len(length(sampleNames)),
                                    pData[[sampling.variable]], maxnum.samples)
        }
        ## update
        mData = mData[,subcols]
        sampleNames = sampleNames[subcols]
        pData = pData[,subcols]
    }


    #feature.selection.quantile = 0.999
    ## subset rows if needed
    if (!is.null(feature.selection.quantile)){
        tmp <- apply(mData,
                     1,
                     get(feature.selection.method),
                     na.rm=TRUE)
        subrows = which(tmp>=quantile(tmp, feature.selection.quantile))
        if (length(subrows)<2){
            stop(sprintf('----- too few rows remained after filtering for top %s most variable rows',
                         feature.selection.quantile))
        }
        #mset.lumi = mset.lumi[subrows,]
        mData = mData[subrows,]
    }

    mData = exprs(mData)

    ## group colors for hclust and/or mds plot
    if ('mds' %in% plots || 'hclust' %in% plots){
        mdist = distFun(distance.method)(t(mData))
        grp_colors=NULL
        if (!is.null(cluster.variable)){
            if (length(cluster.variable)>1){
                cluster.variable = cluster.variable[1]
                warning(sprintf('***** only the first one of [%s] is selected as the group name',
                                paste(cluster.variable, collapse=',')))
            }
            grp_levels = as.character(sort(unique(pData[[cluster.variable]])))
            grp_vector = as.character(pData[[cluster.variable]])
            grp_colors = brewer.pal(n = 8, name = "Dark2")
            if (length(grp_levels)>8){
                grp_colors = rep(GroupColors, ceiling(length(grp_levels)/8))
            }
            grp_colors = grp_colors[seq_len(length(grp_levels))]
            names(grp_colors) = grp_levels
            #col_vector = grp_colors[grp_vector]
        }

        if(!is.null(grp_colors)){
            tmpcols = grp_colors[grp_vector]
        } else {
            tmpcols = 'black'
        }

        ## highlight samples e.g. outliers
        shapes = rep(1, ncol(mData))
        names(shapes) = colnames(mData)
        if (!is.null(highlights)){
            if (all(highlights %in% colnames(mData))){
                shapes[highlights] = 16
            } else {
                warning(sprintf('----- samples/assays in highlights are not a subset of sample names of input data',
                                paste(highlights, sep=',')))
            }
        }

        if ('mds' %in% plots){
            pdf(sprintf ('%s/QC03%s_MDS.pdf',
                         outdir,
                         ifelse(is.null(outname),
                                '',
                                paste('_', outname, sep=''))),
                width=5,
                height=5)
            #print (cmds)
            cmds = cmdscale(mdist, k=2)
            list2return$cmds = cmds
            plot(cmds[,1],
                 cmds[,2],
                 xlab="Dimension 1",
                 ylab="Dimension 2",
                 col=tmpcols,
                 pch=shapes,
                 lwd=1.5,
                 main='MDS of mData')
            if(!is.null(grp_colors)){
                legend('top', grp_levels, text.col=grp_colors, bty='o', cex=1.0, ncol=length(grp_levels)) #horiz=T)
            }
            dev.off()
        }
        #colorCluster = function(hc, labs=NULL, lab.cols=NULL, tip.shapes=NULL, show.labs=TRUE, main.txt=''){
        if ('hclust' %in% plots){
            pdf(sprintf ('%s/QC03%s_hclust.pdf',
                         outdir,
                         ifelse(is.null(outname),
                                '',
                                paste('_', outname, sep=''))),
                width=5,
                height=5)
            hclust.obj = hclust(mdist, method =cluster.method)
            list2return$hclust = hclust.obj #$merge
            colorCluster(hclust.obj,
                         labs=sampleNames,
                         lab.cols=tmpcols,
                         tip.shapes=shapes,
                         show.labs=TRUE,
                         main.txt='clustering of mData')
            if(!is.null(grp_colors)){
                legend('top',
                       grp_levels,
                       text.col=grp_colors,
                       bty='o',
                       bg='white',
                       cex=1.0,
                       ncol=length(grp_levels)) #horiz=T)
            }
            dev.off()
        }
    }

    ### for heatmap3
    if ('heatmap' %in% plots){
        if (!is.null(sidebar.variables)){
            # set colors
            ColSideColors = NULL
            variableColors = list()
            #GroupLevelColorList = NULL # list
            for ( variable in sidebar.variables){
                grp_levels = sort(unique(as.character(pData[[variable]])))
                grp_vector = as.character(pData[[variable]])
                ## colors
                grp_colors = brewer.pal(n = 8, name = "Dark2")
                if (length(grp_levels)>8){
                    grp_colors = rep(grp_colors, ceiling(length(grp_levels)/8))
                }
                grp_colors = grp_colors[seq_len(length(grp_levels))]
                names(grp_colors) = grp_levels
                variableColors[[variable]] = grp_colors
                #for ( mycol in brewColors){
                #    GroupLevelColorList[[variable]][[brewColors[mycol]]] = mycol
                #}
                col_vector = grp_colors[grp_vector] #rep('dark grey', nrow(pData(mset.lumi)))
                #for (g in grp_vector){
                #    col_vector[which(pData(mset.lumi)[[variable]]==g)] = grp_colors[g]
                #}
                ColSideColors = cbind(ColSideColors, col_vector)
                colnames(ColSideColors)[ncol(ColSideColors)] = variable
                #rownames(ColSideColors)=NULL
            }
        }

        pdf(sprintf ('%s/QC03%s_heatmap.pdf',
                     outdir, ifelse(is.null(outname),
                                    '',
                                    paste('_', outname, sep=''))),
            width=5,
            height=5)
        if (is.null(ColSideColors)){
            heatmap3(mData,
                     showRowDendro=FALSE,
                     scale=scaled,
                     distfun=distFun(distance.method),
                     method=cluster.method,
                     balanceColor=balanceColor,
                     cexCol=cexCol,cexRow=cexRow,
                     margins=c(col.margin,row.margin),
                     labRow=labRow,
                     labCol=labCol)
        } else {
            grpnames=NULL; colornames=NULL
            for (variable in names(variableColors)){
                grpnames=c(grpnames, sprintf('(%s)', variable)); colornames=c(colornames, 'black')
                grpnames=c(grpnames, names(variableColors[[variable]])); colornames=c(colornames,variableColors[[variable]])
            }
            heatmap3(mData,
                     ColSideColors=ColSideColors,
                     showRowDendro=FALSE,
                     scale=scaled,
                     distfun=distFun(distance.method),
                     method=cluster.method,
                     balanceColor=balanceColor,
                     cexCol=cexCol,
                     cexRow=cexRow,
                     margins=c(col.margin,row.margin),
                     labRow=labRow,
                     labCol=labCol)
            legend(0,
                   1,
                   xpd=TRUE,
                   grpnames,
                   text.col=colornames,
                   bty='n',
                   cex=1.0,
                   ncol=1,
                   adj=c(1,1))
        }
        dev.off()

        #if (cluster.info && !'hclust' %in% plots){ ## current function won't keep the info
        #    list2return$hclust = hclust(distFun(distance.method)(t(mData)), method =cluster.method)
        #}
    }
    return (list2return)
}

##
# I don't think this is used either...
# vioplot = function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL,
#                     horizontal = FALSE, col = NULL, border = NULL, lty = 1,
#                     lwd = 1, rectCol = 'white', colMed = 'red', pchMed = 19,
#                     at, add = FALSE, wex = 1, drawRect = TRUE)
# {
#     if (is(tmp,'list')){
#         datas = x } else {
#             datas <- list(x, ...)
#         }
#     n <- length(datas)
#     ###
#     for (s in c('col', 'border', 'rectCol', 'colMed')){
#         if (length(get(s))<n){
#             assign(s, rep(get(s), ceiling(n/length(get(s))))[seq_len(n)])
#             #print(s)
#             #print(get(s))
#         }
#     }
#     #if (length(border)<n){
#     #    border = rep(border, ceiling(length(border)/n))[1:n]
#     #}
#     ###
#     if (missing(at))
#         at <- 1:n
#     upper <- vector(mode = "numeric", length = n)
#     lower <- vector(mode = "numeric", length = n)
#     q1 <- vector(mode = "numeric", length = n)
#     q3 <- vector(mode = "numeric", length = n)
#     med <- vector(mode = "numeric", length = n)
#     base <- vector(mode = "list", length = n)
#     height <- vector(mode = "list", length = n)
#     baserange <- c(Inf, -Inf)
#     args <- list(display = "none")
#     if (!(is.null(h)))
#         args <- c(args, h = h)
#     for (i in seq_len(n)) {
#         data <- datas[[i]]
#         data.min <- min(data)
#         data.max <- max(data)
#         q1[i] <- quantile(data, 0.25)
#         q3[i] <- quantile(data, 0.75)
#         med[i] <- median(data)
#         iqd <- q3[i] - q1[i]
#         upper[i] <- min(q3[i] + range * iqd, data.max)
#         lower[i] <- max(q1[i] - range * iqd, data.min)
#         est.xlim <- c(min(lower[i], data.min), max(upper[i],
#                                                    data.max))
#         smout <- do.call("sm.density", c(list(data, xlim = est.xlim),
#                                          args))
#         hscale <- 0.4/max(smout$estimate) * wex
#         base[[i]] <- smout$eval.points
#         height[[i]] <- smout$estimate * hscale
#         t <- range(base[[i]])
#         baserange[1] <- min(baserange[1], t[1])
#         baserange[2] <- max(baserange[2], t[2])
#     }
#     if (!add) {
#         xlim <- if (n == 1)
#             at + c(-0.5, 0.5)
#         else range(at) + min(diff(at))/2 * c(-1, 1)
#         if (is.null(ylim)) {
#             ylim <- baserange
#         }
#     }
#     if (is.null(names)) {
#         label <- seq_len(n)
#     }
#     else {
#         label <- names
#     }
#     boxwidth <- 0.05 * wex
#     if (!add)
#         plot.new()
#     if (!horizontal) {
#         if (!add) {
#             plot.window(xlim = xlim, ylim = ylim)
#             axis(2)
#             axis(1, at = at, label = label)
#         }
#         box()
#         for (i in seq_len(n)) {
#             polygon(c(at[i] - height[[i]],
#                       rev(at[i] + height[[i]])),
#                     c(base[[i]], rev(base[[i]])),
#                     col = col[i],
#                     border = border[i],
#                     lty = lty, lwd = lwd)
#             if (drawRect) {
#                 lines(at[c(i, i)],
#                       c(lower[i], upper[i]),
#                       lwd = lwd,
#                       lty = lty)
#                 rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2,
#                      q3[i], col = rectCol[i])
#                 points(at[i], med[i], pch = pchMed, col = colMed[i])
#             }
#         }
#     }
#     else {
#         if (!add) {
#             plot.window(xlim = ylim, ylim = xlim)
#             axis(1)
#             axis(2, at = at, label = label)
#         }
#         box()
#         for (i in seq_len(n)) {
#             polygon(c(base[[i]], rev(base[[i]])),
#                     c(at[i] - height[[i]], rev(at[i] + height[[i]])),
#                     col = col[i], border = border[i],
#                     lty = lty, lwd = lwd)
#             if (drawRect) {
#                 lines(c(lower[i], upper[i]),
#                       at[c(i, i)],
#                       lwd = lwd,
#                       lty = lty)
#                 rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] +
#                          boxwidth/2, col = rectCol[i])
#                 points(med[i], at[i], pch = pchMed, col = colMed[i])
#             }
#         }
#     }
#     invisible(list(upper = upper,
#                    lower = lower,
#                    median = med,
#                    q1 = q1,
#                    q3 = q3))
# }

## data is a row of a matrix/data.frame, and group is a list matching the data for class (e.g. disease/control),
## type of which is either a factor or a vector
# @importFrom grDevices colorRampPalette
# @importFrom RColorBrewer brewer.pal
## I don't think this is actually used anywhere...
# beeViolin = function(data, group, file=NULL, main='', xlab='', ylab='B-values',
#                      group.colors=NULL){
#     stopifnot(identical(length(data), length(group)))
#     stopifnot(is.numeric(data))
#     ngrps = length(unique(group))
#     if (is.factor(group)){
#         group.labels = levels(group)
#     }else{
#         group.labels = sort(unique(group))
#     }
#
#     if (is.null(group.colors)){
#         group.colors = colorRampPalette(brewer.pal(max(3,ngrps),
#                                                    'Dark2'))(ngrps)
#     } else if (length(group.colors) < ngrps) {
#         group.colors = rep(group.colors,
#                            ceiling(ngrps/length(group.colors)))[seq_len(ngrps)]
#     } else if (length(group.colors) > ngrps) {
#         group.colors = group.colors[seq_len(ngrps)]
#     }
#
#     names(group.labels) = group.colors
#     mycolors = names(group.labels[group])
#
#     outfile = ifelse(is.null(file),
#                      sprintf('violin_beeswarm_plot_%s.pdf', probe),
#                      file)
#
#     ext = substr(outfile,
#                  regexpr("\\.([[:alnum:]]+)$", outfile)+1,
#                  nchar(outfile))
#     if (!ext %in% c('pdf', 'jpeg', 'tiff', 'png')){
#         ext = 'pdf'
#         outfile = paste(outfile, '.pdf', sep='')
#     }
#
#     get(ext)(outfile)
#     beeswarm(data~group,
#              pch=1,
#              pwcol=mycolors,
#              xlab=xlab,
#              ylab=ylab,
#              labels=group.labels,
#              main=main,
#              cex=0.5)
#     vioplot(lapply(group.labels,
#                    function(x) data[group==x]),
#             col=NULL,
#             rectCol='lightblue',
#             border=group.colors,
#             colMed='red',
#             add=TRUE) #
#     dev.off()
# }

##
sampleIndexByPercentageGroup<-function(grpVector, percentage=0.6) {
    elements = unique(grpVector)
    stopifnot(length(elements)==2)
    idx.0<-which(grpVector==elements[1])
    idx.1<-which(grpVector==elements[2])
    size.0<-round(length(idx.0)*percentage)
    size.1<-round(length(idx.1)*percentage)
    sample0.idx<-sample(x=idx.0,size=size.0)
    sample1.idx<-sample(x=idx.1,size=size.1)
    return(c(sample0.idx,sample1.idx))
}

