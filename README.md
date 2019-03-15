## MethyLiution

In order to facilitate meaningful comparisons across the quickly growing number of microarray-based methylation experiments, we have developed an ordered set of procedures which we deploy as a pipeline for the rapid, consistent and user-friendly evaluation of data quality for common formats of these experiments. This pipeline also supports the capability to detect and implement corrective adjustments for notable common inconsistencies arising from experimental and human error.

### Required data

Data should be prepared for analysis by simply moving all related .IDAT files into a single directory.

In addition, a sample sheet/metadata file should be supplied. This must contain at least 4 columns, in the folowing order:

1. Assay ID (derived from standard IDAT filename, minus colour channel label)
2. Sample ID (an arbitrary ID assigned by the experimenter)
3. Experimental group
4. Sex (either F, M, or U if unknown/unavailable)

## Typical operation

The pipeline is divided into a series of 10 steps, which are run in order:

1. Parse data
2. Filter SNP probes and imprinting genes
3. Check that metadata for samples matches inferences
4. Choose the best to carry forwards to subsequent analyses
5. Detect (but do not remove) outliers
6. Colour balancing for 2-colour microarrays
7. Apply inter-assay normalisation
8. Apply Beta-Mixture Quantile normalisation
9. Remove assays with missing metadata
10. Apply surrogate variable analysis to automatically find hidden variables

In most cases, we recommend to simply run the analysis procedure in one pass, like so:

```
input_dir = "C:/path/to/IDAT_files"
sample_sheet = paste(input_dir, "meta.csv", sep = "/")
output_dir = "C:/path/for/output_files"
runPipeline(datadir = input_dir, 
            array = "450K",
            metafile = sample_sheet, 
            outdir = output_dir)
```

By default, we do not run step 10 (surrogate variable analysis) unless explicitly requested, on the grounds that this step is likely to require a more flexible and context-specific analysis than can be easily written into a pipeline function. 

### Included datasets

The list of valid probes (for both EPIC and 450K arrays) required for step 2 is included as an accessible dataset, `probe_data`. This dataset was derived from the list of probes available in the BioConductor packages `IlluminaHumanMethylationEPICanno.ilm10b2.hg19` and `IlluminaHumanMethylation450kprobe`, with filtering applied to remove probes located within 10bp of a SNP in order to avoid misleading results arising from polymorphism-induced hybridisation issues. 

The GEO series [GSE86831](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86831), an EPIC microarray dataset, has been used to produce the plots and outputs presented in the vignette. 
