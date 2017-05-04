#
# Code to compute a fitness metric for a given
#
require(bigWig)

run.cutmask <- function(cutmask, seqoutbias.args, prefix = "run_", bw.split = FALSE, cleanup = TRUE, sqcmd = "seqOutBias") {
    if (!grepl("C", cutmask)) {
        # add C to middle of cutmask
        n = nchar(cutmask)
        stopifnot(n %% 2 == 0)
        left = substr(cutmask, 1, n / 2)
        right = substr(cutmask, n / 2 + 1, n)
        cutmask = paste(left, right, sep='C')

        print(cutmask)
    }
    # build system command
    #
    #   seqOutBias <fasta-file> <bam-file>... [options]
    #
    #  --kmer-mask=<cutmask>
    #  --skip-bed
    #  if bw.split => --stranded
    #  --out=<... filename 1 ...>
    #  --bw=<... filename 2 ...>
    #
    outfilename = paste(prefix, cutmask, ".tbl", sep='')
    out_arg = paste("--out=", outfilename, sep='')
    bwfilename = paste(prefix, cutmask, ".bw", sep='')
    bw_arg = paste("--bw=", bwfilename, sep='')
    mask_arg = paste("--kmer-mask=", cutmask, sep='')
    cmd = paste(sqcmd, seqoutbias.args, mask_arg, "--skip-bed", out_arg, bw_arg)

    # execute
    print(cmd)
    system(cmd)
    

    # clean up - remove files specific to this execution
    if (cleanup) {
        cat("deleting stuff ....\n")
        cat("rm", outfilename,"\n")
        unlink(outfilename)
    }

    # return output file names
    if (bw.split) {
        c(paste(prefix, cutmask, "_plus.bw", sep=''), paste(prefix, cutmask, "_minus.bw", sep=''))
    } else {
        c(bwfilename, bwfilename)
    }
}

#' Evaluate cutmask metric
#'
#' @param motif.sets List of data.frames containing BED6 format coordinates for motif locations; one data.frame per reference transcription factor
#' @param bw.plus BigWig file containing the cutmask scaled data to use for plus strand motif sites
#' @param bw.minus BigWig file containing the cutmask scaled data to use for minus strand motif sites
#'
#' @return metric value
eval.cutmask <- function(motif.sets, bw.plus, bw.minus) {
    metric = 0
    for (bed6 in motif.sets) {
        x = colSums(bed6.step.bpQuery.bigWig(bw.plus, bw.minus, bed6, step = 1, 
                              with.attributes = FALSE, as.matrix = TRUE, follow.strand = TRUE))
        metric = metric + sqrt(sum((x - mean(x))^2))
    }

    return(metric)
}
