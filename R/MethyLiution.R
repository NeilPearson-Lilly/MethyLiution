# library(zeallot)

## meta
assayColNum   <- 1
sampleColNum  <- 2
groupColNum   <- 3
sexColNum     <- 4
## probe
probeColNum   <- 1

# Document any bundled datasets here

#' Filtered probe lists for EPIC and 450K datasets
#'
#' The list of valid probes (for both EPIC and 450K arrays) required for step 2 is included as an accessible dataset, `probe_data`. This dataset was derived from the list of probes available in the BioConductor packages `IlluminaHumanMethylationEPICanno.ilm10b2.hg19` and `IlluminaHumanMethylation450kprobe`, with filtering applied to remove probes located within 10bp of a SNP identified in dbSNP build 140 in order to avoid misleading results arising from polymorphism-induced hybridisation issues.
#'
#' @docType data
#' @keywords datasets
#' @name probe_data
#' @usage data(probe_data)
#' @format A list of two data frames.
NULL



#' Set column names
#'
#' Set column names of input dataset to consistent values.
#' @param mset.lumi A methylation dataset object
#' @return New column names
#' @importFrom Biobase pData pData<-
#' @importFrom AnnotationDbi colnames
#' @importFrom logging loginfo

setColumnNames <- function(mset.lumi) {
    for (s in c('assayColNum', 'sampleColNum', 'sexColNum', 'groupColNum')){
        assign(sub('Num', 'Name', s), colnames(pData(mset.lumi))[get(s)])
        loginfo(sprintf('----- %s is [%s]', sub('Num', 'Name', s),
                        colnames(pData(mset.lumi))[get(s)]))
    }
    if(length(unique(pData(mset.lumi)[[sampleColName]]))==1){
        pData(mset.lumi)[[sampleColName]] = paste(seq_len(nrow(pData(mset.lumi))),
                                                  pData(mset.lumi)[[sampleColName]],
                                                  sep=':')
    }
    col.name <- list(assayColName=assayColName,
                     sampleColName=sampleColName,
                     sexColName=sexColName,
                     groupColName=groupColName)
    col.name
}

#' Step 1: Read dataset files
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param datadir Path to directory containing IDAT files
#' @param metafile Path to file containing experimental metadata
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param array Microarray probeset used: valid options are "450K" and "EPIC" (default "450K")
#' @param threads Number of threads to run certain tasks in parallel
#' @param FDB Methylation annotation set (auto-filled by runPipeline() function)
#' @return List - column names and methylumi dataset
#' @importFrom wateRmelon readEPIC
#' @importFrom minfi getAnnotation
#' @importFrom logging loginfo logwarn
#' @importFrom methylumi getBarcodes
#' @importFrom BiocGenerics annotation
#' @importFrom Biobase fData pData pData<- sampleNames
#' @importFrom lumi sampleNames<-
#' @import parallel
#' @import IlluminaHumanMethylationEPICmanifest

readData <- function(datadir, metafile, logfile, mset_file, array = '450K',
                     threads = 1, FDB) {
  loginfo(sprintf('1.1 ##### get the file names from %s', datadir))

  filenames <- list.files(datadir,
                          pattern='*\\.idat',
                          recursive = TRUE,
                          include.dirs = TRUE)
  if(length(filenames)==0){stop("No idat files found!")}
  barcodes <- getBarcodes(filenames)
  nlibs <- length(barcodes)
  if(2*nlibs != length(filenames)){stop("idat files are not complete!")}

  ##barcodes = barcodes[1:20]
  # ???threads=1 ## it won't work in cluster

  loginfo(sprintf('1.2 ##### read raw data from %s files of %s samples',
                  length(filenames),
                  nlibs))
  if (threads>1) {
    loginfo(sprintf('----- parallel in %s threads to read the idat files',
                    threads))
    options(mc.cores=threads)
    mset.lumi <- wateRmelon::readEPIC(oob=TRUE,
                          idatPath=datadir,
                          parallel=TRUE,
                          force=TRUE)
  } else {
    mset.lumi <- wateRmelon::readEPIC(oob=TRUE,
                          idatPath=datadir,
                          parallel=FALSE,
                          force=TRUE)
  }

  if(! grepl(array, annotation(mset.lumi), ignore.case=TRUE)){
    stop("Wrong array type was specified!")
  }

  mset.lumi<-as(mset.lumi, "MethyLumiM")

  if(array == "450K"){
    if (!all(c('COLOR_CHANNEL', 'CHROMOSOME', 'POSITION', 'CN_TYPE') %in%
             colnames(fData(mset.lumi)))){
      loginfo(sprintf('1.3 ##### add probe information %s to mset.lumi',
                      paste(c('"COLOR_CHANNEL',
                              'CHROMOSOME',
                              'POSITION',
                              'CN_TYPE"'),
                            collapse=' ')))
      # Function from utilities.r
      mset.lumi <- addAnnotation(as(mset.lumi, "MethyLumiM"),
                                 lib=FDB,
                                 annotationColumn=c('COLOR_CHANNEL',
                                                    'CHROMOSOME',
                                                    'POSITION',
                                                    'CN_TYPE'))
    }
  }

  if (ncol(pData(mset.lumi)) < 4){
    loginfo(sprintf('1.4 ##### read meta file:: %s', metafile))

    if(!file.exists(metafile)){stop("No meta file found!")}
    meta_dat <- read.csv(metafile, header=TRUE,
                         sep=',',
                         stringsAsFactors=FALSE,
                         na.strings = c("NA" , "." ))
    sampleNames(mset.lumi) <- gsub("^(*.\\/)+", "", sampleNames(mset.lumi))

    if(sum(meta_dat[[assayColNum]] %in% sampleNames(mset.lumi))==0){
      stop("The first column in the meta file is not barcode or the meta file does not include the methylation samples!")
    }

    matches <- match(sampleNames(mset.lumi), meta_dat[[assayColNum]])
    # use assay (data-wise) instead of sample (specimen-wise)

    if (sum(is.na(matches)) > 0){
      logwarn(sprintf ('----- there are **%s** samples missing in meta data',
                       sum(is.na(matches))))
      meta_dat <- do.call(rbind, lapply(seq_len(length(matches)),
                                        function(x) {
                                            if (is.na(matches[x])) {
                                                rep(NA, ncol(meta_dat))}
                                            else {
                                                meta_dat[matches[x],]
                                                }
                                            }))
      meta_dat[[assayColNum]] <- sampleNames(mset.lumi)
    } else {
      meta_dat <- meta_dat[matches,]  #matches[!is.na(matches)]
    }

    #write.csv(meta_dat, metafile_updated, quote=F, row.names=T)
    #write.table(as.character(sampleNames(mset.lumi)), 'samples_lumi.csv', quote=F, row.names=F)

    loginfo(sprintf('1.5 ##### add meta data to mset.lumi as phenoData'))
    rownames(meta_dat) <- sampleNames(mset.lumi)
    pData(mset.lumi) <- meta_dat
  }

  loginfo(sprintf('1.6 ##### get the column names in meta data'))
  # for (s in c('assayColNum', 'sampleColNum', 'sexColNum', 'groupColNum')){
  #   assign(sub('Num', 'Name', s), colnames(pData(mset.lumi))[get(s)])
  #   loginfo(sprintf('----- %s is [%s]', sub('Num', 'Name', s), colnames(pData(mset.lumi))[get(s)]))
  # }
  # if(length(unique(pData(mset.lumi)[[sampleColName]]))==1){
  #   pData(mset.lumi)[[sampleColName]] = paste(1:nrow(pData(mset.lumi)), pData(mset.lumi)[[sampleColName]], sep=':')
  # }
  col.name = setColumnNames(mset.lumi)

  loginfo(sprintf('----- column names are: sample %s; sex %s; group %s',
                  col.name$sampleColName,
                  col.name$sexColName,
                  col.name$groupColName))
  loginfo(sprintf('----- summary information of original mset.lumi'))
  sink(logfile, append=TRUE); show(mset.lumi); sink()

  # col.name <- list(assayColName=assayColName, sampleColName=sampleColName,
  #                  sexColName=sexColName, groupColName=groupColName)

  loginfo(sprintf('1.7 ##### save mset.lumi and col.name'))
  save(mset.lumi, file=gsub('mset.lumi$','QC01_mset.lumi',
                            mset_file))
  save(col.name, file=paste0(mset_file,
                             ".col.name"))
  list(col.name = col.name,
       mset.lumi = mset.lumi)
}

#' Step 2: Filter SNP probes and imprinting genes
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param begin Step the pipeline begins at
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param qc_probe_plot_file Path to file where probe QC information will be written (auto-filled by runPipeline() function)
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param array Microarray probeset used: valid options are "450K" and "EPIC" (default "450K")
#' @param col.name List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)
#' @param mset.lumi Microarray dataset object (auto-filled by runPipeline() function)
#' @return mset.lumi dataset, with SNP probes and imprinting genes filtered out
#' @importFrom logging loginfo
#' @importFrom Biobase fData fData<- pData sampleNames
#' @importFrom grDevices pdf dev.off
#' @importFrom methylumi qc.probe.plot
#' @importFrom ggplot2 theme element_blank element_text labs
#' @importFrom utils data
#' @import IlluminaHumanMethylationEPICanno.ilm10b2.hg19

removeSnpProbesAndImprintingGenes <- function(begin, mset_file,
                                              qc_probe_plot_file, logfile,
                                              array, col.name = NULL,
                                              mset.lumi = NULL) {
    if(begin == 2){
        # These get written at the end of the step 1 function.
        load(gsub('mset.lumi$', 'QC01_mset.lumi', mset_file))
        load(paste0(mset_file, ".col.name"))
    }
    assayColName  <- col.name$assayColName
    sexColName    <- col.name$sexColName
    sampleColName <- col.name$sampleColName
    groupColName  <- col.name$groupColName

    loginfo(sprintf('2.1 ##### remove probes associated with SNPs or imprinting genes according to pre-prepared %s probe set',
    array))

    # probefile = paste0(array, "_PROBES.RData")
    # if(!file.exists(probefile)){ stop("No probe file found!") }
    # ??? where to put probe file?
    # probe_dat <- read.csv(probefile, header=TRUE, stringsAsFactors=FALSE)
    # probe_dat = readRDS(probefile)
    # load(probefile)
    # The BioConductor requirements leave us with no good alternatives to doing
    # this.
    # The data is picked up automatically from sysdata.rda
    # probe_dat = get(paste0("probes_", array))

    # We don't need to do that any more. We now just load a list of dataframes
    # (probe_data) which has the names we need in data frames inside it.
    load("data/probe_data.Rdata")
    # This gives us this nice sanity check back too:
    if(!paste0("probes_", array) %in% names(probe_data)) {
        stop("No probe data for this assay type found!")
    }
    # This is where the filtering happens: all we do is just pick out every ID
    # that's in the given probe list. We don't use the metadata at all.
    # mset.lumi <- mset.lumi[fData(mset.lumi)[,1] %in% probe_dat[[probeColNum]],]
    mset.lumi <- mset.lumi[Biobase::fData(mset.lumi)[,1] %in% probe_data[[paste0("probes_", array)]][['Name']],]

    if(array == "EPIC"){
        # Add annotation
        # This now comes from the IlluminaHumanMethylationEPICanno.ilm10b2.hg19
        # package.
        data(Locations)
        anno_epic = probe_data[[paste0("probes_", array)]]
        anno_epic$probes_EPIC.COLOR_CHANNEL = as.character(anno_epic$probes_EPIC.COLOR_CHANNEL)
        anno_epic = cbind(anno_epic, Locations[c("chr", "pos")][match(anno_epic$Name, rownames(Locations)),])
        names(anno_epic) = c("Name", "COLOR_CHANNEL", "CN_TYPE", "CHROMOSOME",
                             "POSITION")
        anno_epic<-anno_epic[,c("Name", "CHROMOSOME", "POSITION",
                                "COLOR_CHANNEL", "CN_TYPE")]

        anno_epic$COLOR_CHANNEL[anno_epic$COLOR_CHANNEL==""] <- "Both"
        rownames(anno_epic) <- anno_epic$Name
        anno_epic<-anno_epic[rownames(Biobase::fData(mset.lumi)),]
        names(anno_epic) <-c ("Probe_ID", "CHROMOSOME", "POSITION",
                              "COLOR_CHANNEL", "CN_TYPE")
        anno_epic$DESIGN = lapply(anno_epic$COLOR_CHANNEL, function (x) ifelse(x == "Both", "II", "I"))

        # Try to merge the values in anno_epic with the fData in mset.lumi,
        # preserving order. That's easy. Just use rownames from that data frame.
        # Chances are it's the same already but we want to be sure.

        # We just can't use the fData<- setter - it's completely broken and
        # simply will not work in this case. We have to do it the hard way.

        new_fData = fData(mset.lumi)
        new_fData$CHROMOSOME = anno_epic[rownames(new_fData),]$CHROMOSOME
        new_fData$POSITION = anno_epic[rownames(new_fData),]$POSITION
        new_fData$CN_TYPE = anno_epic[rownames(new_fData),]$CN_TYPE

        # fData(mset.lumi) <- anno_epic
        # Now this has to be dealt with by putting the data in the right place
        attr(attr(mset.lumi, 'featureData'), 'data') <- new_fData

        feature_metadata = attr(attr(mset.lumi, 'featureData'), 'varMetadata')
        feature_metadata <- rbind(feature_metadata, "CHROMOSOME" = c("Chromosome"))
        feature_metadata <- rbind(feature_metadata, "POSITION" = c("Location index on chromosome"))
        feature_metadata <- rbind(feature_metadata, "CN_TYPE" = c("Modification type"))
        attr(attr(mset.lumi, 'featureData'), 'varMetadata') <- feature_metadata
    }

    loginfo(sprintf('----- summary information of filtered mset.lumi according to pre-prepared %s probe set',
                    array))
    sink(logfile, append=TRUE); show(mset.lumi); sink()

    loginfo(sprintf('2.2 ##### QC plot:: using control probes'))

    pdf(qc_probe_plot_file)
    print(qc.probe.plot(mset.lumi, by.type=FALSE) +
    theme(axis.text.y=element_blank(),
    axis.title=element_text(face="bold",
                            size="11",
                            color="black"),
                            legend.position="right") +
                            labs(title='Negative & normalisation control probes (raw mset)'))
    dev.off()

    loginfo(sprintf('2.3 ##### save mset.lumi'))
    save(mset.lumi, file=gsub('mset.lumi$',
                              'QC02_mset.lumi',
                              mset_file))
    mset.lumi
}

#' Step 3: Check that metadata for samples matches inferences
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param begin Step the pipeline begins at
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param gender_stringent Should the pipeline discard samples where inferred gender is inconsistent with metadata?
#' @param sex_label_cluster_file Path to file where probe QC information will be written (auto-filled by runPipeline() function)
#' @param outdir Path to output directory
#' @param col.name List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)
#' @param mset.lumi Microarray dataset object (auto-filled by runPipeline() function)
#' @return mset.lumi dataset
#' @importFrom logging loginfo
#' @importFrom Biobase fData pData pData<- sampleNames
#' @importFrom stats cutree
#' @importFrom utils write.csv

checkReportedSex <- function(begin, mset_file, logfile, gender_stringent,
                             sex_label_cluster_file, outdir,
                             col.name = NULL, mset.lumi = NULL) {
    if(begin == 3){
        load(gsub('mset.lumi$', 'QC02_mset.lumi', mset_file))
        load(paste0(mset_file, ".col.name"))
    }
    assayColName  <- col.name$assayColName
    sexColName    <- col.name$sexColName
    sampleColName <- col.name$sampleColName
    groupColName  <- col.name$groupColName

    loginfo(sprintf('3.1 ##### format the sex labels as F/M/U'))
    pData(mset.lumi)[[sexColName]] <- formatSexNames(pData(mset.lumi)[[sexColName]])

    if(length(unique(pData(mset.lumi)[[sexColNum]]))>1){
        loginfo(sprintf('3.2 ##### QC plot:: chrX probes'))
        chrom = 'chrX'
        groupByChrom <- plotHeatmap(mset.lumi[fData(mset.lumi)$CHROMOSOME==chrom,],
                                    outname=chrom,
                                    plots=c('heatmap',
                                            'mds',
                                            'hclust'),
                                    outdir=outdir,
                                    maxnum.samples=NULL,
                                    sampling.variable=NULL,
                                    sidebar.variables=sexColName,
                                    cluster.variable=sexColName,
                                    distance.method="euclidean",
                                    feature.selection.method = 'mad',
                                    feature.selection.quantile = 0.9,
                                    scaled='row',
                                    cluster.method="ward.D2",
                                    balanceColor=TRUE,
                                    cexCol=1, cexRow=1,
                                    col.margin=5, row.margin=5,
                                    labRow=NA, labCol=NA)
        use_hclust <- cutree(groupByChrom$hclust, 2)
        sexes <- as.character(pData(mset.lumi)[[sexColName]])
        sexes_by_cluster <- as.character(reLabelVector(sexes, use_hclust))

        if (!identical(sexes_by_cluster, sexes)){
            logwarn(sprintf('----- %s assay(s) found with in-concordant sex clustering; results are output to "%s"',
                            sum(sexes!=sexes_by_cluster),
                            sex_label_cluster_file))
            write.csv(data.frame(Assay=sampleNames(mset.lumi),
                                 SexLabel=sexes,
                                 SexCluster=sexes_by_cluster,
                                 Warn=ifelse(sexes==sexes_by_cluster, 0, 1)),
                      sex_label_cluster_file,
                      quote=FALSE,
                      row.names=FALSE)
            if (gender_stringent){
                samples <- sampleNames(mset.lumi)[sexes==sexes_by_cluster]
                mset.lumi <- mset.lumi[,!(pData(mset.lumi)[[sampleColName]] %in%
                                              samples)]
                logwarn(sprintf('----- %s assay(s) with suspicious sex labelling being removed',
                                sum(sexes!=sexes_by_cluster)))
            }
        }

    } else {
        loginfo(sprintf('----- Only one sex found, no clustering by sex performed.'))
    }

    logwarn(sprintf('3.3 ##### remove the probes of chrX and Y'))
    mset.lumi = mset.lumi[!(fData(mset.lumi)$CHROMOSOME %in% c("chrX","chrY")),]

    loginfo(sprintf('----- summary information of filtered mset.lumi to remove chrX/Y probes%s',
                    ifelse(gender_stringent,
                           ' and potentially mis-labelled samples',
                           '')))
    sink(logfile, append=TRUE); show(mset.lumi); sink()

    loginfo(sprintf('3.4 ##### save mset.lumi'))
    save(mset.lumi, file=gsub('mset.lumi$',
                              'QC03_mset.lumi',
                              mset_file))
    mset.lumi
}

#' Step 4: Where multiple arrays are present for some samples, choose the best to carry forwards to subsequent analyses
#'         Detect the duplicated assays for the same sample, and keep the one assay with the lowest mean detection p-value
#'         Remove methylation values with higherdetection p-value (larger than --maxDetectionPval)
#'         Remove probes with failed detection (p-value> maxDetectionPval) in more than -maxDetectionProp number of samples
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param begin Step the pipeline begins at
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param max_pval Upper threshold for considering a result significant
#' @param max_prop For each probe, the freq/proportion of assays that fail by max_pval of detection
#' @param detection_freq_infofile Duplicated assays output file (auto-filled by runPipeline() function)
#' @param detection_filter_infofile Detection failure table by probe file (auto-filled by runPipeline() function)
#' @param outdir Path to output directory
#' @param col.name List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)
#' @param mset.lumi Microarray dataset object (auto-filled by runPipeline() function)
#' @return mset.lumi dataset, with less performant assays for each sample removed
#' @importFrom logging loginfo
#' @importFrom Biobase fData pData sampleNames
#' @importFrom lumi beta2m
#' @importFrom utils write.csv
#' @importFrom methylumi exprs<- detection betas

chooseBestArray <- function(begin, mset_file, logfile, max_pval, max_prop,
                            detection_freq_infofile, detection_filter_infofile,
                            outdir, col.name = NULL, mset.lumi = NULL) {
    # mset.lumi, assayColName, sampleColName
    if(begin == 4){
        load(gsub('mset.lumi$',
                  'QC03_mset.lumi',
                  mset_file))
        load(paste0(mset_file, ".col.name"))
    }
    assayColName  <- col.name$assayColName
    sexColName    <- col.name$sexColName
    sampleColName <- col.name$sampleColName
    groupColName  <- col.name$groupColName

    assays_to_del = NULL
    probes_to_del = NULL

    flt_df <- data.frame(sampleColName=pData(mset.lumi)[[sampleColName]],
                         assayColName=pData(mset.lumi)[[assayColName]])
    # names(flt_df) <- sapply(names(flt_df), get)  # This is nice but seems to generate more trouble than it's worth.
    names(flt_df) = c(sampleColName, assayColName)

    if (sum(duplicated(pData(mset.lumi)[[sampleColName]], incomparables=NA))){
        loginfo(sprintf(paste('4.1 ##### mark the duplicated assays to remove',
                              ' - keep the one with the least mean detection ',
                              'p-val across all remained probes')))
        dup_samples = unique(pData(mset.lumi)[[sampleColName]][duplicated(pData(mset.lumi)[[sampleColName]],
                                                                          incomparables=NA)])
        for (dup in dup_samples){
            dupAssays <- as.character(pData(mset.lumi)[[assayColName]][which(pData(mset.lumi)[[sampleColName]] == dup)])
            avePvals <- sort(apply(detection(mset.lumi)[,dupAssays], 2, mean))
            # named vect
            assays_to_del <- c(assays_to_del,
                               names(avePvals[-1]))
        }

        flt_df[['DeDuplicated']] <- ifelse(flt_df[[assayColName]] %in%
                                               assays_to_del, 'YES', 'NO')

        loginfo(sprintf('----- %s duplicated assay(s) of %s sample(s) being removed; see file %s for details',
                        length(assays_to_del),
                        length(assays_to_del),
                        detection_filter_infofile,
                        outdir))
    }

    loginfo(sprintf('4.2 ##### remove assay-probes with detection p-val>%s',
                    max_pval))

    filter_by_detection = is.na(exprs(mset.lumi)[which(detection(mset.lumi) >
                                                           max_pval,
                                                       arr.ind=TRUE)]) = TRUE

    loginfo(sprintf('----- %s calls are removed as bad detection; see files "%s" and "%s" for details',
                    sum(is.na(exprs(mset.lumi))),
                    detection_filter_infofile,
                    detection_freq_infofile))

    loginfo(sprintf('----- output table of number of probes (including duplicates if any) with absent calls for each assay to %s',
                    detection_filter_infofile))
    flt_df[['NumAbsentProbes']] = apply(beta2m(betas(mset.lumi)),
                                        2,
                                        function(x) sum(is.na(x)))
    write.csv(flt_df, detection_filter_infofile, quote=FALSE, row.names=FALSE)

    assays_to_del <- unique(assays_to_del)
    if (!is.null(assays_to_del)){
        loginfo(sprintf('4.3 ##### remove %s duplicated assays identified according to detection pvals from mset.lumi',
                        length(assays_to_del)))
        if(nrow(pData(mset.lumi)) == length(assays_to_del)){
            stop("All arrays were removed!")
        }

        mset.lumi = mset.lumi[,!(pData(mset.lumi)[[assayColName]] %in%
                                     assays_to_del)]
        loginfo(sprintf('----- summary information of filtered mset.lumi after removing duplicates'))
        sink(logfile, append=TRUE); show(mset.lumi); sink()
    }

    loginfo(sprintf('4.4 #### calculate frequency of failed detection in all remained assays for each probe'))
    nas_by_probe <- apply(lumi::beta2m(methylumi::betas(mset.lumi)),
                          1,
                          function(x) sum(is.na(x)))
    probes_to_del <- c(probes_to_del,
                       which(nas_by_probe >= max_prop*length(sampleNames(mset.lumi))))

    loginfo(sprintf('4.5 #### output table of number of remained assays with absent calls for each probe to %s',
                    detection_freq_infofile))
    table_na_by_probe <- data.frame(table(nas_by_probe))
    names(table_na_by_probe)=c('AssayFrequency', 'NumAbsentProbes')
    table_na_by_probe[['HiFreqFailure']] = 'NO'
    table_na_by_probe[['HiFreqFailure']][as.numeric(as.character(table_na_by_probe[['AssayFrequency']]))>=
                                             max_prop*length(sampleNames(mset.lumi))] = 'YES'

    write.csv(table_na_by_probe,
              detection_freq_infofile,
              quote=FALSE,
              row.names=FALSE)

    if ((!is.null(probes_to_del)) & length(probes_to_del)>0) {
        loginfo(sprintf('4.6 ##### remove %s probes with high frequency of detection failures from mset.lumi',
                        length(probes_to_del)))
        mset.lumi = mset.lumi[-probes_to_del,]
        loginfo(sprintf(paste('----- summary information of filtered mset.lumi',
                               'after removing high-frequency-failure probes')))
        sink(logfile, append=TRUE); show(mset.lumi); sink()
    }

    loginfo(sprintf('4.7 ##### save mset.lumi'))
    save(mset.lumi, file=gsub('mset.lumi$',
                              'QC04_mset.lumi',
                              mset_file))
    mset.lumi
}

#' Step 5: Detect (but do not remove) outliers
#' Makes a hierarchical clustering plot of all samples.
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param begin Step the pipeline begins at
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param outliers_file Path to outliers output file
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param max_plot Maximum number of plots to create at this step
#' @param outdir Path to output directory
#' @param col.name List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)
#' @param mset.lumi Microarray dataset object (auto-filled by runPipeline() function)
#' @return list - mset.lumi dataset, breakdown of sub-assays
#' @importFrom logging loginfo
#' @importFrom utils write.csv

detectOutliers <- function(begin, mset_file, outliers_file, logfile,
                           max_plot, outdir,
                           col.name = NULL, mset.lumi = NULL) {
    # mset.lumi, groupColName
    if(begin == 5){
        load(gsub('mset.lumi$',
                  'QC04_mset.lumi',
                  mset_file))
        load(paste0(mset_file, ".col.name"))
    }
    assayColName  <- col.name$assayColName
    sexColName    <- col.name$sexColName
    sampleColName <- col.name$sampleColName
    groupColName  <- col.name$groupColName

    loginfo(sprintf('5.1 ##### QC plots and outlier detection'))

    max_plot <- min(ncol(mset.lumi), max_plot)
    clustering <- plotQC(mset.lumi,
                         outdir=outdir,
                         step='05',
                         scaled=TRUE,
                         plots=c('color',
                                 'density',
                                 'outlier'),
                         maxnum.samples=max_plot,
                         sampling.variable=groupColName,
                         cluster.variable=groupColName,
                         distance.method='correlation',
                         cluster.method="average",
                         outlier.threshold=2,
                         outlier.plot=TRUE)

    sub_assays <- clustering[['subset.samples']]

    if(is.null(clustering$dist.outliers)){
        clustering$dist.outliers<-data.frame(distThreshold=rep(FALSE,
                                                               ncol(mset.lumi)),
                                             distCluster=rep(FALSE,
                                                             ncol(mset.lumi)),
                                             outlier=rep(FALSE,
                                                         ncol(mset.lumi)))}

    write.csv(data.frame(Assay=sampleNames(mset.lumi),
                         clustering$dist.outliers),
              outliers_file, quote=FALSE, row.names=FALSE)
    loginfo(sprintf('----- %s potential outlier(s) detected; if any, see %s for details',
                    sum(clustering$dist.outliers$outlier),
                    outliers_file))

    loginfo(sprintf('5.2 ##### save mset.lumi'))
    # save(sub_assays, file=sub_assays_file)
    save(sub_assays, file=paste0(mset_file, ".sub.assays"))
    list(mset.lumi = mset.lumi,
         sub_assays = sub_assays)
}

#' Step 6: Colour balancing for 2-colour microarrays
#' Applies a normalisation method to remove systemic differences between colour channels.
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param begin Step the pipeline begins at
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param outdir Path to output directory
#' @param col.name List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)
#' @param mset.lumi Microarray dataset object (auto-filled by runPipeline() function)
#' @param sub_assays Object containing breakdown of sub-assays, made during Step 5.
#' @return list - mset.lumi dataset, colour-balanced dataset object
#' @importFrom logging loginfo
#' @importFrom lumi lumiMethyC

colourAdjustment <- function(begin, mset_file, logfile, outdir,
                             col.name = NULL, mset.lumi = NULL,
                             sub_assays = NULL) {
    if(begin == 6){
        load(gsub('mset.lumi$',
                  'QC04_mset.lumi',
                  mset_file))
        load(paste0(mset_file, ".col.name"))
        # load(paste0(mset_file, ".sub.assays"))
        # load(paste0(sub_assays_file, ".sub_assays"))
    }
    assayColName  <- col.name$assayColName
    sexColName    <- col.name$sexColName
    sampleColName <- col.name$sampleColName
    groupColName  <- col.name$groupColName

    loginfo(sprintf('6.1 ##### color adjustment'))
    lumi.C <- lumiMethyC(mset.lumi) #method='quantile' quantile color balance adj
    if (is.null(sub_assays)){
        plotQC_C <- plotQC(lumi.C,
                           outdir=outdir,
                           step='06',
                           outname='afterColorAdj',
                           scaled=TRUE,
                           plots=c('color', 'density'))
    } else {
        loginfo(sprintf('----- %s samples/assays selected for plotting purpose',
                        length(sub_assays)))
        plotQC_C <- plotQC(lumi.C[, sub_assays],
                           outdir=outdir,
                           step='06',
                           outname='afterColorAdj',
                           scaled=TRUE,
                           plots=c('color', 'density'))
    }

    loginfo(sprintf('6.2 ##### save mset.lumi and lumi.C'))
    save(lumi.C,    file=gsub('mset.lumi$',
                              'QC06_mset.lumi.C',
                              mset_file))
    list(mset.lumi = mset.lumi, lumi.C = lumi.C)
}

#' Step 7: Apply inter-assay normalisation
#' Apply a broader inter-array normalisation procedure than the colour balancing in Step 6, determined by user input.
#' Options are quantile, SSN, and none.
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param begin Step the pipeline begins at
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param normalisation_method Normalisation method: one of 'quantile', 'ssn', or 'none'.
#' @param outdir Path to output directory
#' @param col.name List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)
#' @param mset.lumi Microarray dataset object (auto-filled by runPipeline() function)
#' @param lumi.C Colour-balanced dataset object (auto-filled by runPipeline() function)
#' @param sub_assays Object containing breakdown of sub-assays, made during Step 5. (auto-filled by runPipeline() function)
#' @return list - mset.lumi dataset, meth.B normalised dataset object (exprs), meth.M normalised dataset object (betas)
#' @importFrom logging loginfo
#' @importFrom lumi lumiMethyN
#' @importFrom Biobase exprs
#' @importFrom methylumi betas

normalisation <- function(begin, mset_file, logfile, normalisation_method,
                          outdir,
                          col.name = NULL, mset.lumi = NULL, lumi.C = NULL,
                          sub_assays = NULL) {
    if(begin == 7){
        load(gsub('mset.lumi$', 'QC06_mset.lumi.C', mset_file))
        load(paste0(mset_file, ".col.name"))
        # load(paste0(mset_file, ".sub.assays"))
    }
    assayColName  <- col.name$assayColName
    sexColName    <- col.name$sexColName
    sampleColName <- col.name$sampleColName
    groupColName  <- col.name$groupColName

    if (normalisation_method %in% c('quantile', 'ssn')){
        loginfo(sprintf('7.1 ##### normalise two-color combined raw data using %s method',
                        normalisation_method))
        lumi.N <- lumiMethyN(lumi.C, method=normalisation_method)
        meth.M = exprs(lumi.N)
        meth.B = betas(lumi.N)
        if (is.null(sub_assays)){
            plotQC_N = plotQC(lumi.N,
                              outdir=outdir,
                              step='07',
                              outname='afterNorm',
                              scaled=TRUE,
                              plots=c('color', 'density'))
        } else {
            plotQC_N = plotQC(lumi.N[, sub_assays],
                              outdir=outdir,
                              step='07',
                              outname='afterNorm',
                              scaled=TRUE,
                              plots=c('color', 'density'))
        }
    } else {
        logwarn(sprintf('7.1 ##### no normalisation - %s',
                        normalisation_method))
        meth.M = exprs(lumi.C)
        meth.B = betas(lumi.C)
    }

    loginfo(sprintf('7.2 ##### save mset.lumi, meth.M and meth.B'))
    save(meth.B,
         file=gsub('mset.lumi$', 'QC07_mset.lumi.meth.B',
                   mset_file))
    save(meth.M,
         file=gsub('mset.lumi$', 'QC07_mset.lumi.meth.M',
                   mset_file))
    list(mset.lumi = mset.lumi,
         meth.B = meth.B,
         meth.M = meth.M)
}

#' Step 8: Apply Beta-Mixture Quantile normalisation, producing BMIQ objects
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param begin Step the pipeline begins at
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param BMIQ_dir Path to a directory where BMIQ normalised data can be written to files (auto-filled by runPipeline() function)
#' @param outdir Path to output directory
#' @param threads Number of threads to run certain tasks in parallel
#' @param col.name List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)
#' @param mset.lumi Microarray dataset object (auto-filled by runPipeline() function)
#' @param meth.B Object containing normalised data, made during Step 7. (auto-filled by runPipeline() function)
#' @return list - mset.lumi dataset, BMIQ-normalised data object (B values), BMIQ-normalised data object (M-values)
#' @importFrom logging loginfo
#' @importFrom doParallel registerDoParallel
#' @importFrom wateRmelon BMIQ
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%

quantileNormalisation <- function(begin, mset_file, logfile, BMIQ_dir,
                                  outdir, threads,
                                  col.name = NULL, mset.lumi = NULL,
                                  meth.B = NULL) {
    if(begin == 8){
        load(gsub('mset.lumi$',
                  'QC04_mset.lumi',
                  mset_file)) # no need
        load(gsub('mset.lumi$',
                  'QC07_mset.lumi.meth.B',
                  mset_file))
        load(paste0(mset_file, ".col.name"))
    }
    assayColName  <- col.name$assayColName
    sexColName    <- col.name$sexColName
    sampleColName <- col.name$sampleColName
    groupColName  <- col.name$groupColName

    probe_type_vector = ifelse(fData(mset.lumi)[['COLOR_CHANNEL']] == 'Both',
                               2,
                               1)

    loginfo(sprintf('8.1 ##### BMIQ normalisation'))
    workdir = getwd()
    dir.create(BMIQ_dir, showWarnings = FALSE)
    setwd(BMIQ_dir)
    loginfo(sprintf('----- has changed directory to %s for BMIQ plots',
                    getwd()))

    bmiq.B = matrix(NA,
                    nrow=nrow(meth.B),
                    ncol=ncol(meth.B))

    cl <- makeCluster(max(min(threads, detectCores()-1, 20), 1))
    # what is the optimal number

    # suppressPackageStartupMessages(require(doParallel))
    # suppressPackageStartupMessages(require(RPMM))
    registerDoParallel(cl, cores = cl)

    bmiq.B <- foreach(i = seq_len(ncol(meth.B)),
                      .packages = c("wateRmelon"),
                      .combine = cbind) %dopar% {
        try({
            tmpBMIQ <- BMIQ(meth.B[,i],
                            probe_type_vector,
                            nL = 3,
                            doH = TRUE,
                            nfit = 50000,
                            th1.v = c(0.1, 0.76),
                            th2.v = c(0.17, 0.7),
                            niter= 5,
                            tol = 0.001,
                            plots = TRUE,
                            sampleID = colnames(meth.B)[i],
                            pri=FALSE)
            Sys.sleep(3)
            tmpDf <- data.frame(tmpBMIQ$nbeta)
            colnames(tmpDf) = colnames(meth.B)[i]
            tmpDf

        })
    }

    stopifnot(identical(colnames(bmiq.B),
                        colnames(meth.B)))
    stopCluster(cl)

    rownames(bmiq.B) <- rownames(meth.B)

    setwd(workdir)
    loginfo(sprintf('----- done with BMIQ; change back to previous working directory %s',
                    getwd()))
    #outputProbeLevelData('bmiq.B', outdir, 'QC08') # no need

    loginfo(sprintf('8.2 ##### convert B-values to M-values'))
    bmiq.M <- convertBM(bmiq.B)
    #outputProbeLevelData('bmiq.M', outdir, 'QC08') # no need

    loginfo(sprintf('8.3 ##### save mset.lumi and bmiq.B'))
    save(bmiq.B,    file=gsub('mset.lumi$', 'QC08_mset.lumi.bmiq.B',
                              mset_file))
    save(bmiq.M,    file=gsub('mset.lumi$', 'QC08_mset.lumi.bmiq.M',
                              mset_file))
    list(mset.lumi = mset.lumi,
         bmiq.B = bmiq.B,
         bmiq.M = bmiq.M)
}

#' Step 9: Remove assays with missing metadata
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param begin Step the pipeline begins at
#' @param metafile Path to file containing experimental metadata
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param feature_data_file Path to a file containing feature data (auto-filled by runPipeline() function)
#' @param missing_meta_file Path to a file where assays with missing data are explicitly listed (auto-filled by runPipeline() function)
#' @param meta_data_file Path to a file containing valid metadata from samples that have passed this filter (auto-filled by runPipeline() function)
#' @param outdir Path to output directory
#' @param col.name List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)
#' @param mset.lumi Microarray dataset object (auto-filled by runPipeline() function)
#' @param meth.B Object containing normalised data, made during Step 7. (auto-filled by runPipeline() function)
#' @param meth.M Object containing normalised data, made during Step 7. (auto-filled by runPipeline() function)
#' @param bmiq.B Object containing BMIQ-normalised data, made during Step 8. (auto-filled by runPipeline() function)
#' @param bmiq.M Object containing BMIQ-normalised data, made during Step 8. (auto-filled by runPipeline() function)
#' @return list - mset.lumi dataset, BMIQ-normalised data object (B values)
#' @importFrom logging loginfo
#' @importFrom utils write.csv read.csv

removeAssaysWithMissingMetadata <- function(begin, mset_file, logfile,
                                            feature_data_file, metafile,
                                            missing_meta_file, meta_data_file,
                                            outdir, col.name = NULL,
                                            mset.lumi = NULL,
                                            meth.B = NULL, meth.M = NULL,
                                            bmiq.B = NULL, bmiq.M = NULL) {
    if(begin == 9){
        load(gsub('mset.lumi$',
                  'QC04_mset.lumi',
                  mset_file))
        load(gsub('mset.lumi$',
                  'QC07_mset.lumi.meth.B',
                  mset_file))
        load(gsub('mset.lumi$',
                  'QC07_mset.lumi.meth.M',
                  mset_file))
        load(gsub('mset.lumi$',
                  'QC08_mset.lumi.bmiq.B',
                  mset_file))
        load(gsub('mset.lumi$',
                  'QC08_mset.lumi.bmiq.M',
                  mset_file))
        load(paste0(mset_file, ".col.name"))
    }
    assayColName  <- col.name$assayColName
    sexColName    <- col.name$sexColName
    sampleColName <- col.name$sampleColName
    groupColName  <- col.name$groupColName

    loginfo(sprintf('9.1 ##### remove assays with missing data regarding %s',
                    groupColName))
    assays_missing_meta <- union(which(is.na(pData(mset.lumi)[[groupColName]])),
                                 which(pData(mset.lumi)[[groupColName]] %in%
                                           c('NA', '.', '', 'na')))

    loginfo(sprintf('----- output the assay names with missing meta data to %s',
                    missing_meta_file))
    write.csv(data.frame(assaysNoMeta=sampleNames(mset.lumi)[assays_missing_meta]),
              missing_meta_file,
              quote=FALSE,
              row.names=FALSE)

    loginfo(sprintf('9.2 ##### update the meth.B, meth.M, mset.lumi and bmiq.M'))
    if(length(assays_missing_meta)>0){
        meth.M <- meth.M[, -assays_missing_meta]
        meth.B <- meth.B[, -assays_missing_meta]
        bmiq.M <- bmiq.M[, -assays_missing_meta]
        bmiq.B <- bmiq.B[, -assays_missing_meta]
        mset.lumi <- mset.lumi[, -assays_missing_meta]
    }

    ##
    #' @importFrom utils write.csv
    outputProbeLevelData = function(dat_name_string, outdir, outname,
                                    samples = NULL,
                                    probes = NULL){
        result_file = sprintf('%s/%s%s.csv',
                              outdir,
                              ifelse(is.null(outname),
                                     '',
                                     paste(outname, '_', sep='')),
                              dat_name_string)
        dat = get(dat_name_string)
        if (is.null(probes)) {
            probes = rownames(dat)
        }
        if (is.null(samples)){
            samples = colnames(dat)
        }
        write.csv(data.frame(probeID = probes,
                             dat[probes, samples],
                             check.names = FALSE),
                  file = result_file,
                  row.names = FALSE,
                  quote = FALSE)
        loginfo(sprintf('----- output %s to %s',
                        dat_name_string, result_file))
    }

    loginfo(sprintf('9.3 ##### output meth.B, meth.M, bmiq.B and bmiq.M to files'))
    probe_ids <- as.character(rownames(meth.B))
    outputProbeLevelData('meth.B', outdir, 'QC09')
    outputProbeLevelData('meth.M', outdir, 'QC09')
    outputProbeLevelData('bmiq.B', outdir, 'QC09')
    outputProbeLevelData('bmiq.M', outdir, 'QC09')

    loginfo(sprintf('----- output array feature data to %s',
                    feature_data_file))
    write.csv(fData(mset.lumi),
              file = feature_data_file,
              row.names = FALSE,
              quote = FALSE)


    loginfo(sprintf('----- output updated meta data to %s', meta_data_file))
    meta_dat <- read.csv(metafile,
                         header = TRUE,
                         sep= ',',
                         stringsAsFactors = FALSE,
                         na.strings = c("NA" , "." ))
    matches  <- match(sampleNames(mset.lumi), meta_dat[[assayColNum]])
    meta_dat <- meta_dat[matches,]
    write.csv(meta_dat,
              file = meta_data_file,
              row.names = FALSE,
              quote = FALSE)


    loginfo(sprintf('9.4 ##### save mset.lumi and bmiq.M'))
    save(mset.lumi,
         file=gsub('mset.lumi$',
                   'QC09_mset.lumi',
                   mset_file))
    save(bmiq.M,
         file=gsub('mset.lumi$',
                   'QC09_mset.lumi.bmiq.M',
                   mset_file))
    list(mset.lumi = mset.lumi,
         bmiq.M = bmiq.M)
}

#' Step 10: Apply surrogate variable analysis to automatically find hidden variables in the data
#'
#' This probably isn't something you want to run directly. Just use the runPipeline() function.
#'
#' @keywords internal
#' @param begin Step the pipeline begins at
#' @param mset_file Path to file in which dataset object is stored between steps (auto-filled by runPipeline() function)
#' @param sv_file Path to file where any surrogate variables found are written as output (auto-filled by runPipeline() function)
#' @param logfile Path to log file (auto-filled by runPipeline() function)
#' @param col.name List of column names, generated by setColumnNames() (auto-filled by runPipeline() function)
#' @param mset.lumi Microarray dataset object (auto-filled by runPipeline() function)
#' @param bmiq.M Object containing BMIQ-normalised data, made during Step 8. (auto-filled by runPipeline() function)
#' @return mset.lumi dataset
#' @importFrom logging loginfo
#' @importFrom sva num.sv sva
#' @importFrom utils write.csv
#' @importFrom stats model.matrix

findSurrogateVariables <- function(begin, mset_file, sv_file, logfile,
                                   col.name = NULL, mset.lumi = NULL,
                                   bmiq.M = NULL) {
    if(begin == 10){
        load(gsub('mset.lumi$',
                  'QC09_mset.lumi',
                  mset_file))
        load(gsub('mset.lumi$',
                  'QC09_mset.lumi.bmiq.M',
                  mset_file))
        load(paste0(mset_file, ".col.name"))
    }
    assayColName <- col.name$assayColName
    sexColName <- col.name$sexColName
    sampleColName <- col.name$sampleColName
    groupColName <- col.name$groupColName

    loginfo(sprintf('10 ##### SVA using updated mset.lumi'))

    mod <- model.matrix(~as.factor(pData(mset.lumi)[[groupColName]]))
    mod0 <- model.matrix(~1, data=pData(mset.lumi))
    loginfo(sprintf('----- estimate n.sv'))
    n.sv <- num.sv(bmiq.M,
                   mod,
                   method='be')
    loginfo(sprintf('----- n.sv is %s', n.sv))

    loginfo(sprintf('----- sva'))
    svobj <- sva(as.matrix(bmiq.M ),
                 mod,
                 mod0,
                 n.sv = n.sv)
    loginfo(sprintf('----- output sv to a table file %s',
                    sv_file))
    write.csv(data.frame(Sample=sampleNames(mset.lumi),
                         svobj$sv),
              sv_file,
              quote=FALSE,
              row.names=FALSE)

    loginfo(sprintf('##### Done with probe-level data processing!'))
    mset.lumi
}

#' runPipeline
#'
#' This function is the means of accessing quality control functions of the MethyLiution package. Most procedures are invoked in a set order due to requirements for output from previous steps, so a single function invoking those procedures as a fixed workflow is preferable. \cr
#' Analysis is broken down into a pipeline of 10 steps: \cr
#' \itemize{
#'  \item 1. Read data. \cr
#'  \item 2. Remove SNP probes and imprinting genes. \cr
#'  \item 3. Check reported sex against inferred sex. \cr
#'  \item 4. Remove all but the most performant assay per sample. \cr
#'  \item 5. Detect (but do not remove) outliers. \cr
#'  \item 6. Normalise for dye colour. \cr
#'  \item 7. Apply inter-assay normalisation. \cr
#'  \item 8. Apply beta-mixture quantile (intra-sample) normalisation. \cr
#'  \item 9. Remove assays with missing data. \cr
#'  \item 10. Apply surrogate variable analysis to discover hidden variables. \cr
#' }
#' For more details, see vignette(MethyLiution).
#'
#' @param datadir Path to directory containing IDAT files
#' @param metafile Path to file containing experimental metadata
#' @param outdir Path to output directory
#' @param array Microarray probeset used: valid options are "450K" and "EPIC" (default "450K")
#' @param max_pval Upper threshold for considering a result significant (default 0.05)
#' @param max_prop For each probe, the freq/proportion of assays that fail by max_pval of detection (default 0.5)
#' @param max_plot Maximum number of plot output files to create (default 30)
#' @param threads Number of threads to run certain tasks in parallel (default 1)
#' @param gender_stringent Should the pipeline discard samples where inferred gender is inconsistent with metadata? (default FALSE)
#' @param normalisation_method Normalisation method: one of 'quantile', 'ssn', or 'none'. (default 'quantile')
#' @param begin Step the pipeline begins at (default 1)
#' @param end Step the pipeline ends at (default 9)
#' @return path to output directory
#' @importFrom logging loginfo basicConfig addHandler writeToFile
#' @importFrom parallel detectCores
#' @importFrom zeallot %<-%
#' @examples runPipeline()
#' @export

runPipeline <- function(datadir, metafile, outdir, array = '450K',
                        max_pval = 0.05, max_prop = 0.5, max_plot = 30,
                        threads = 1, gender_stringent = FALSE,
                        normalisation_method = 'quantile',
                        begin = 1, end = 9) {
    if (threads > 1) {
        ncores = max(min(detectCores()-1, 20), 1)
        threads = ifelse(threads>ncores,
                         ncores,
                         threads)
    }

    ## user defined parameters
    # probefile        <- file.path(array, 'PROBES.csv')

    if(begin > end) {
        stop("--end should be equal or larger than --begin!")
    }
    if(max_pval > 1 | max_pval<0) {
        stop("--maxDetectionPval should be within 0:1")
    }
    if(max_prop > 1 | max_prop<0) {
        stop("--maxDetectionProp should be within 0:1")
    }
    if(max_plot <= 0) {
        max_plot <- NULL
    }

    if(array=="450K") {
        FDB="FDb.InfiniumMethylation.hg19"
    } else if(array=="EPIC") {
        FDB="IlluminaHumanMethylationEPICanno.ilm10b2.hg19"
    } else {
        stop("Only HM450 and EPIC arrays are supported!")
    }

    ## for test only
    #datadir = 'rawdata'
    #outdir = 'PLEVEL1'
    #probefile = 'srcfiles//PROBES.csv'
    #metafile = 'srcfiles//TEST_META.csv'
    #max_pval = 0.05
    #max_plot = 30
    #threads =6
    #gender_stringent = FALSE
    #normalisation.method = 'median'
    #begin = 4
    #end = 4
    ##

    dir.create(outdir, showWarnings = FALSE)

    #### define file names
    # logfile                   <- sprintf('%s/%s/Run_%s.log', getwd(), outdir,
    #                                      format(Sys.time(),
    #                                             '%y%m%d%H%M')) #need full path!
    logfile                   <- sprintf('%s/Run_%s.log',
                                         normalizePath(outdir),
                                         format(Sys.time(),
                                                '%y%m%d%H%M')) #need full path!
    mset_file                 <- sprintf('%s/%s',
                                         outdir,
                                         'mset.lumi')
    qc_probe_plot_file        <- sprintf('%s/QC02_plot_control_probes.pdf',
                                         outdir)
    sex_label_cluster_file    <- sprintf('%s/QC03_Sex_Clusters.csv',
                                         outdir)
    detection_filter_infofile <- sprintf('%s/QC04_Duplicated_Assays.csv',
                                         outdir)
    detection_freq_infofile   <- sprintf('%s/QC04_Detection_Failure_Table_By_Probe.csv',
                                         outdir)
    outliers_file             <- sprintf('%s/QC05_Outliers.csv',
                                         outdir)
    BMIQ_dir                  <- sprintf('%s/QC08BMIQ_PLOTS',
                                         outdir)
    missing_meta_file         <- sprintf('%s/QC09_assays_missing_meta.csv',
                                         outdir)
    feature_data_file         <- sprintf('%s/QC09_Array_Features.csv',
                                         outdir)
    sv_file                   <- sprintf('%s/QC10_SVs.csv',
                                         outdir)
    meta_data_file            <- sprintf('%s/meta_updated.csv',
                                         outdir)
    mset.lumi = NULL
    lumi.C = NULL
    meth.B = NULL
    meth.M = NULL

    basicConfig(level='DEBUG')
    addHandler(writeToFile,
               file=logfile,
               level='DEBUG')

    if(begin > 1){
        load(paste0(mset_file,
                    ".col.name"))
    }
    if(begin > 5){
        load(paste0(mset_file,
                    ".sub.assays"))
    }

    if (begin <= 1 & end >= 1){
        c(col.name, mset.lumi) %<-% readData(datadir = datadir,
                                             metafile = metafile,
                                             logfile = logfile,
                                             mset_file = mset_file,
                                             array = array,
                                             threads = threads,
                                             FDB = FDB)
    }
    if (begin <= 2 & end >= 2){
        mset.lumi = removeSnpProbesAndImprintingGenes(begin = begin,
                                                      mset_file = mset_file,
                                                      qc_probe_plot_file = qc_probe_plot_file,
                                                      logfile = logfile,
                                                      array = array,
                                                      col.name = col.name,
                                                      mset.lumi = mset.lumi)
    }
    if (begin <= 3 & end >= 3){
        mset.lumi = checkReportedSex(begin = begin,
                                     mset_file = mset_file,
                                     logfile = logfile,
                                     gender_stringent = gender_stringent,
                                     sex_label_cluster_file = sex_label_cluster_file,
                                     outdir = outdir,
                                     col.name = col.name,
                                     mset.lumi = mset.lumi)
    }
    if (begin <= 4 & end >= 4){
        mset.lumi = chooseBestArray(begin = begin,
                                    mset_file = mset_file,
                                    logfile = logfile,
                                    max_pval = max_pval,
                                    max_prop = max_prop,
                                    detection_freq_infofile = detection_freq_infofile,
                                    detection_filter_infofile = detection_filter_infofile,
                                    outdir = outdir,
                                    col.name = col.name,
                                    mset.lumi = mset.lumi)
    }
    if (begin <= 5 & end >= 5){
        c(mset.lumi, sub_assays) %<-% detectOutliers(begin = begin,
                                                     mset_file = mset_file,
                                                     outliers_file = outliers_file,
                                                     logfile = logfile,
                                                     max_plot = max_plot,
                                                     outdir = outdir,
                                                     col.name = col.name,
                                                     mset.lumi = mset.lumi)
    }
    if (begin <= 6 & end >= 6){
        c(mset.lumi, lumi.C) %<-% colourAdjustment(begin = begin,
                                                   mset_file = mset_file,
                                                   logfile = logfile,
                                                   outdir = outdir,
                                                   col.name = col.name,
                                                   mset.lumi = mset.lumi,
                                                   sub_assays = sub_assays)
    }
    if (begin <= 7 & end >= 7){
        c(mset.lumi, meth.B, meth.M) %<-%  normalisation(begin = begin,
                                                         mset_file = mset_file,
                                                         logfile = logfile,
                                                         normalisation_method = normalisation_method,
                                                         outdir = outdir,
                                                         col.name = col.name,
                                                         mset.lumi = mset.lumi,
                                                         lumi.C = lumi.C,
                                                         sub_assays = sub_assays)
    }
    if (begin <= 8 & end >= 8){
        c(mset.lumi, bmiq.B, bmiq.M) %<-% quantileNormalisation(begin = begin,
                                                                mset_file = mset_file,
                                                                logfile = logfile,
                                                                BMIQ_dir = BMIQ_dir,
                                                                outdir = outdir,
                                                                threads = threads,
                                                                col.name = col.name,
                                                                mset.lumi = mset.lumi,
                                                                meth.B = meth.B)
    }
    if (begin <= 9 & end >= 9){
        c(mset.lumi, bmiq.M) %<-% removeAssaysWithMissingMetadata(begin = begin,
                                                                  mset_file = mset_file,
                                                                  logfile = logfile,
                                                                  metafile = metafile,
                                                                  feature_data_file = feature_data_file,
                                                                  missing_meta_file = missing_meta_file,
                                                                  meta_data_file = meta_data_file,
                                                                  outdir = outdir,
                                                                  col.name = col.name,
                                                                  mset.lumi = mset.lumi,
                                                                  meth.B = meth.B,
                                                                  meth.M = meth.M,
                                                                  bmiq.B = bmiq.B,
                                                                  bmiq.M = bmiq.M)
    }
    if (begin <= 10 & end >= 10){
        mset.lumi = findSurrogateVariables(begin = begin,
                                           mset_file = mset_file,
                                           sv_file = sv_file,
                                           logfile = logfile,
                                           col.name = col.name,
                                           mset.lumi = mset.lumi,
                                           bmiq.M = bmiq.M)
    }

    outdir
}


#' Cite MethyLiution
#'
#' If this software has proved useful to you, please cite it!
#'
#' @return Citation for this software, in BibTeX format
#' @examples cite_me()
#' @export

cite_me <- function() {
    c = "No citeable publication yet!"
    c
}
