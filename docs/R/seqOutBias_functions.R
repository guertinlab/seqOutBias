library(lattice)
library(bigWig)
library(gplots) 
library(gtools)

#functions
composites.test.naked <- function(path.dir, composite.input, region=20, step=1, grp = 'DNase') { 
  vec.names = c('chr','start','end')
  hmap.data = list()
  composite.df=data.frame(matrix(ncol = 6, nrow = 0))
  for (mod.bigWig in Sys.glob(file.path(path.dir, "*.bigWig"))) {
      factor.name = strsplit(strsplit(mod.bigWig, "/")[[1]][length(strsplit(mod.bigWig, "/")[[1]])], '\\.')[[1]][1]
      print(factor.name)
      vec.names = c(vec.names, factor.name)
      wiggle = load.bigWig(mod.bigWig)
      bpquery = window.step(composite.input, wiggle, region, step)
      subsample = subsampled.quantiles.metaprofile((bpquery[[1]]))
      #alternative functions: bootstrapped.confinterval.metaprofile, 
      #confinterval.metaprofile, subsampled.quantiles.metaprofile
      mult.row = ncol(bpquery[[1]])
      hmap.data[[factor.name]] = bpquery[[1]]
      df.up <- data.frame(matrix(ncol = 6, nrow = mult.row))
      df.up[, 1] <- colMeans(bpquery[[1]])
      df.up[, 2] <- seq((-1 * region) + 0.5 * step, region - 0.5 * step, by = step)
      df.up[, 3] <- matrix(data = factor.name, nrow=mult.row, ncol=1)
      df.up[, 4] <- subsample$top
      df.up[, 5] <- subsample$bottom
      df.up[, 6] <- matrix(data = grp, nrow=mult.row, ncol=1)
      composite.df = rbind(composite.df, df.up)
      unload.bigWig(wiggle)
  }
  colnames(composite.df) <- c('est', 'x', 'cond', 'upper', 'lower', 'grp')
  composite.df = composite.df[composite.df[,2] >= -1000 & composite.df[,2] <= 1000,]
  for (cond in (1:length(hmap.data))) {
      rownames(hmap.data[[cond]]) = paste(composite.input[,1], ':', composite.input[,2], '-', composite.input[,3], sep='')
      colnames(hmap.data[[cond]]) = seq((-1 * region) + 0.5 * step, region - 0.5 * step, by = step)
  }
  return(list(composite.df, hmap.data))
}

randomRows <- function(df,n){
   return(df[sample(nrow(df),n),])
}

processes.comps.not.hotspots <- function(fimo = fimo.df, path.dir = 'UW_DNase', region = 200, fctr.name = 'efl1') {
    factor.name = fctr.name
    all.fimo.bed = fimo[fimo[,2] > region,]
    chrom.bW = load.bigWig(Sys.glob(file.path(path.dir, "*.bigWig"))[1])
   for (i in 1:length(chrom.bW$chroms)) {
        all.fimo.bed = subset(all.fimo.bed, (all.fimo.bed[,1] == chrom.bW$chroms[i]) & 
                                  (all.fimo.bed[,3] < chrom.bW$chromSizes[i]) | all.fimo.bed[,1] != chrom.bW$chroms[i])         
    }
    unload.bigWig(chrom.bW)
    fimo.composites = composites.test.naked(path.dir, all.fimo.bed, region = region, grp = factor.name)
    return(fimo.composites)
}

window.step <- function(bed, bigWig, halfWindow, step) {
  windowSize = (2*halfWindow) %/% step
  midPoint = floor((as.numeric(as.character(bed[,2])) + as.numeric(as.character(bed[,3]))) / 2)
  start = (midPoint - halfWindow)
  end = start + windowSize*step
  if ((as.numeric(as.character(bed[1,2])) + as.numeric(as.character(bed[1,3]))) %% 2 == 0) {
    bed[,2] = start
    bed[,3] = end
  } else {
    bed[,2] = start
    bed[,3] = end
    bed[,2][bed[,6] == '-'] = bed[,2][bed[,6] == '-'] + 1
    bed[,3][bed[,6] == '-'] = bed[,3][bed[,6] == '-'] + 1
  }
  matrix.comp = bed6.step.bpQuery.bigWig(bigWig, bigWig, bed, 1, op = "avg", follow.strand = TRUE)
  res = do.call(rbind, matrix.comp)
  return(list(res, matrix.comp))
}

bedTools.intersect<-function(functionstring="/usr/local/bin/bedtools/intersectBed",bed1,bed2,opt.string="") {
  a.file=tempfile()
  b.file=tempfile()
  out   =tempfile()
  options(scipen =99)
  write.table(bed1,file=a.file,quote=F,sep="\t",col.names=F,row.names=F)
  write.table(bed2,file=b.file,quote=F,sep="\t",col.names=F,row.names=F)
  command=paste(functionstring, opt.string,"-a",a.file,"-b",b.file,">",out,sep=" ")
  cat(command,"\n")
  try(system(command))
  res=read.table(out,header=F, comment.char='')
  unlink(a.file);unlink(b.file);unlink(out)
  return(res[,c(4,5,6,7,8,9)])
}

my.panel.bands <-
    function(x, y, upper, lower,
             fill, col,
             subscripts, ..., font, fontface)
{
    upper <- upper[subscripts]
    lower <- lower[subscripts]
    panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                  col = fill, border = FALSE,
                  ...)
}

composites.func.panels.naked.chromatin <- function(dat, fact = 'Factor', summit = 'Summit', num=90, 
              col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
              rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), rgb(1/2,1/2,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
  count = length(unique(dat$grp))
  ct.cons = 0
  lst.cons = list()
  unique(dat$grp)[order(unique(dat$grp))]
  for (i in unique(dat$grp)[order(unique(dat$grp))]) {
      ct.cons= ct.cons + 1
      lst.cons[[ct.cons]] = c(min(dat[dat$grp == i,]$lower), max(dat[dat$grp == i,]$upper))
  }
  pdf(paste('composite_', fact, '_signals_', summit, '_peaks.pdf', sep=''), width=3.43, 
      height=ceiling((count)) * 3.00) 
  print(xyplot(est ~ x|grp, group = cond, data = dat,
               type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
               xlim=c(-(num),(num)),
               ylim = lst.cons,
               col = col.lines,
               auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,2,2,2,2,2), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
               cex.axis=1.0,
               par.strip.text=list(cex=0.9, font=1, col='black'),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               #lwd=2,
               ylab = list(label = paste(fact," Cut Frequency", sep=''), cex =0.8),
               xlab = list(label = paste("Distance from ", summit, " center",sep=''), cex =0.8),
               upper = dat$upper,
               fill = fill.poly,
               lower = dat$lower,
               strip = function(..., which.panel, bg) {
                 bg.col = c("grey85")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
             },
               panel = function(x, y, ...){
                   panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                   panel.xyplot(x, y, ...)
       }
  ))
  dev.off()
}

parse.fimo <- function(file) {
  fimo.data = read.table(file, colClasses = c('character','character','integer','integer',
                                   'character','numeric', 'numeric', 'character'))
  res = fimo.data[,c(2,3,4,6,7,5,1)]
  colnames(res) = c('chr', 'start', 'end', 'score', 'pval', 'strand', 'motif')
  return(res)
}


parse.mast <- function(file, motif.num = 1) {
  mast.data = read.table(file, colClasses = c('character','character','integer','integer','numeric','numeric'))
  filename = paste(strsplit(file, '.txt')[[1]][1], '.bed', sep='')
  print(filename)
  chrom = vector(mode="character", length = nrow(mast.data))
  start = vector(mode="integer", length = nrow(mast.data))
  end = vector(mode="integer", length = nrow(mast.data))
  strand.motif = vector(mode="integer", length = nrow(mast.data))
  score = vector(mode="numeric", length = nrow(mast.data))
  pval = vector(mode="numeric", length = nrow(mast.data))

  for (j in 1:nrow(mast.data)) {
    chrom[j] = strsplit(mast.data[j,1], ":")[[1]][1]
    start[j] = as.numeric(strsplit(strsplit(mast.data[j,1], ":")[[1]][2], "-")[[1]][1]) + mast.data[j,3]
    end[j] = as.numeric(strsplit(strsplit(mast.data[j,1], ":")[[1]][2], "-")[[1]][1]) + mast.data[j,4]
    strand.motif[j] = as.character(mast.data[j,2])
    score[j] = mast.data[j,5]
    pval[j] = mast.data[j,6]
  }
  
  arg = data.frame(cbind(chrom, start, end, score, pval, strand.motif))
  arg = arg[arg$strand.motif == paste('+', motif.num, sep='') | arg$strand.motif ==  paste('-', motif.num, sep=''),]

  strand = vector(mode="character", length = nrow(arg))
  motif = vector(mode="character", length = nrow(arg))
  for (j in 1:nrow(arg)) {
    strand[j] = strsplit(as.character(arg[j,6]), "")[[1]][1]
    motif[j] = strsplit(as.character(arg[j,6]), "")[[1]][2]
  }
  res = cbind(arg[,c(T,T,T,T,T,F)], strand, motif)
  colnames(res) = c('chr', 'start', 'end', 'score', 'pval', 'strand', 'motif')
  res$start = as.numeric(as.character(res$start))
  res$end = as.numeric(as.character(res$end))

  write.table(res, file=filename, sep='\t', quote=F, row.names=F, col.names=F)
  return(res)
}

cycle.fimo.new.not.hotspots <- function(path.dir.fimo = FALSE, path.dir.mast = FALSE, rand.rows = 100000,
                                        path.dir.bigWig = '~/UW_DNase_MCF7', window = 200, exp = 'ATAC') {
  composite.df=data.frame(matrix(ncol = 6, nrow = 0))
  #first thing to do is to convert either mast or fimo input to a motif data frame
  if (path.dir.fimo != FALSE & path.dir.mast != FALSE) {
      print('choose either fimo input or mast input, not both')
     }
  if (path.dir.fimo == FALSE & path.dir.mast == FALSE) {
      print('choose either fimo input or mast input')
     }
  if (path.dir.fimo ==FALSE & path.dir.mast != FALSE) {
      for (mast.file in Sys.glob(file.path(path.dir.mast, "*_mast.txt"))) {
          factor.name = strsplit(strsplit(mast.file, "/")[[1]]
              [length(strsplit(mast.file, "/")[[1]])], '\\_mast')[[1]][1]
          print(factor.name)
          parsed = parse.mast(mast.file)
          tf.x = processes.comps.not.hotspots(fimo = parsed,
              path.dir =  path.dir.bigWig, region = window, fctr.name = factor.name)
          composite.df = rbind(composite.df, tf.x[[1]])
      } 
      composites.func.panels.naked.chromatin(composite.df, exp, summit = 'Motif', num = window - (window *0.2))
  }
  if (path.dir.mast == FALSE & path.dir.fimo != FALSE) {
      for (fimo.file in Sys.glob(file.path(path.dir.fimo, "*_fimo.txt"))) {
          factor.name = strsplit(strsplit(fimo.file, "/")[[1]]
              [length(strsplit(fimo.file, "/")[[1]])], '\\_fimo')[[1]][1]
          print(factor.name)
          parsed = parse.fimo(fimo.file)
#fimo output can be hundreds of thousands of motifs and overwhelm R--we take a random subset:          
          parsed = randomRows(parsed[grep('_',parsed[,1], invert =TRUE),], rand.rows)
          tf.x = processes.comps.not.hotspots(fimo = parsed,
              path.dir =  path.dir.bigWig, region = window, fctr.name = factor.name)
          composite.df = rbind(composite.df, tf.x[[1]])
      } 
      composites.func.panels.naked.chromatin(composite.df, exp, summit = 'Motif', num = window - (window *0.2))
  }
  return(composite.df)
}


composites.func <- function(dat, fact = 'Factor', summit = 'Summit', num=90, 
              col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
              rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), rgb(1/2,1/2,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
  count = length(unique(dat$grp))
  ct.cons = 0
  lst.cons = list()
  unique(dat$grp)[order(unique(dat$grp))]
  for (i in unique(dat$grp)[order(unique(dat$grp))]) {
      ct.cons= ct.cons + 1
      lst.cons[[ct.cons]] = c(min(dat[dat$grp == i,]$lower), max(dat[dat$grp == i,]$upper))
  }
  pdf(paste('composite_', fact, '_signals_', summit, '_peaks.pdf', sep=''), width=8.2, 
      height=5.43) 
  print(xyplot(est ~ x|grp, group = cond, data = dat,
               type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
               xlim=c(-(num),(num)),
               ylim = lst.cons,
               col = col.lines,
               auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,0.8,2,0.8,0.8,0.8), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
               cex.axis=1.0,
               par.strip.text=list(cex=0.9, font=1, col='black'),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               layout=c(4,2),
               ylab = list(label = paste(fact," Detection Frequency", sep=''), cex =0.8),
               xlab = list(label = paste("Distance from ", summit, " center",sep=''), cex =0.8),
               upper = dat$upper,
               fill = fill.poly,
               lower = dat$lower,
               strip = function(..., which.panel, bg) {
                 bg.col = c("grey85")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
             },
               panel = function(x, y, ...){
                   panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                   panel.xyplot(x, y, ...)
                   panel.abline(v = 9.5, lty = 2, col = "grey40")
                   panel.abline(v = -9.5, lty = 2, col = "grey40")

       }
  ))
  dev.off()
}

revcomp <- function(x) {
  return(step1 = paste(rev(unlist(strsplit(chartr("ATGC","TACG", toupper(x)), split=""))), collapse=""))
}

plot.barchart.lig <- function(df.barchart, filename = "temp.pdf", w = 12, h = 4) {
    rat = df.barchart[,6]/df.barchart[,5]
    pdf(filename, width=w, height=h)
    print(barchart((rat~df.barchart[,8]),
                   ylim = c(0, max(rat[is.finite(rat)], na.rm =TRUE)+ 0.01 * max(rat[is.finite(rat)], na.rm =TRUE)),
                   col='grey85',
                   ylab=paste('ratio',df.barchart[1,9], 'to alt. seq. orientation', sep = ' '),
                   xlab = 'alternative sequence orientation',
                   scales=list(x=list(rot=45)),
                   panel=function(...) {
                       panel.barchart(...)
                       panel.abline(h=1, lty= 2, col = 'grey40')
                   }
                   ))
    dev.off()
}

matrix.func <- function(counts.table) {
    ligation = cbind(counts.table, substring(counts.table[,2],1,3))
    ligation[,8] = apply(ligation, 1, function(row) revcomp(row[7]))
    ligation[,9] = substring(counts.table[,2],4,6)
    colnames(ligation) = c(colnames(counts.table), 'V7', 'rcKmerUp','KmerDown')
    
    mat = data.frame(matrix(nrow=64, ncol= 64))
    count = 0
    for (mer in unique(ligation$KmerDown)) {
        count = count + 1
        temp = ligation[ligation$KmerDown == mer,]
        rto = temp[,6]/temp[,5]
        mat[,count] = rto
    }
    colnames(mat) = unique(ligation$KmerDown)
    mat = do.call(data.frame,lapply(mat, function(x) replace(x, is.infinite(x),NA)))
    rownames(mat) = unique(temp[,8])
    mat = mat[order(rownames(mat)) , order(colnames(mat))]
    mat = as.matrix(mat)
    return(mat)
}


composites.func.pro <- function(dat, fact = 'Factor', summit = 'Summit', num=90, 
              col.lines = c(rgb(0,0,1,1/2), rgb(1,0,0,1/2),  rgb(0.1,0.5,0.05,1/2), rgb(0,0,0,1/2),  
              rgb(1/2,0,1/2,1/2), rgb(0,1/2,1/2,1/2), rgb(1/2,1/2,0,1/2)), fill.poly = c(rgb(0,0,1,1/4), 
              rgb(1,0,0,1/4), rgb(0.1,0.5,0.05,1/4),rgb(0,0,0,1/4), rgb(1/2,0,1/2,1/4))) {
  count = length(unique(dat$grp))
  ct.cons = 0
  lst.cons = list()
  unique(dat$grp)[order(unique(dat$grp))]
  for (i in unique(dat$grp)[order(unique(dat$grp))]) {
      ct.cons= ct.cons + 1
      lst.cons[[ct.cons]] = c(min(dat[dat$grp == i,]$lower), max(dat[dat$grp == i,]$upper))
  }
  pdf(paste('composite_', fact, '_signals_', summit, '_peaks.pdf', sep=''), width=10.2, 
      height=10.2) 
  print(xyplot(est ~ x|grp, group = cond, data = dat,
               type = 'l',
               scales=list(x=list(cex=0.8,relation = "free"), y =list(cex=0.8, relation="free")),
               xlim=c(-(num),(num)),
               ylim = lst.cons,
               col = col.lines,
               auto.key = list(points=F, lines=T, cex=0.8),
               par.settings = list(superpose.symbol = list(pch = c(16), col=col.lines, cex =0.5), 
                   superpose.line = list(col = col.lines, lwd=c(2,2,0.8,2,0.8,0.8,0.8), 
                       lty = c(1,1,1,1,1,1,1,1,1))),
               cex.axis=1.0,
               par.strip.text=list(cex=0.9, font=1, col='black'),
               aspect=1.0,
               between=list(y=0.5, x=0.5),
               layout=c(3,3),
               ylab = list(label = paste(fact," Cut Frequency", sep=''), cex =0.8),
               xlab = list(label = paste("Distance from ", summit, " center",sep=''), cex =0.8),
               upper = dat$upper,
               fill = fill.poly,
               lower = dat$lower,
               strip = function(..., which.panel, bg) {
                 bg.col = c("grey85")
                 strip.default(..., which.panel = which.panel, bg = rep(bg.col, length = which.panel)[which.panel])
             },
               panel = function(x, y, ...){
                   panel.superpose(x, y, panel.groups = 'my.panel.bands', ...)
                   panel.xyplot(x, y, ...)
       }
  ))
  dev.off()
}

pswm.func <- function(x.ligation, out = 'outfilename') {
    a = lapply(strsplit(as.character(x.ligation), ''), "[", 1)
    b = lapply(strsplit(as.character(x.ligation), ''), "[", 2)
    c = lapply(strsplit(as.character(x.ligation), ''), "[", 3)
    d = lapply(strsplit(as.character(x.ligation), ''), "[", 4)
    e = lapply(strsplit(as.character(x.ligation), ''), "[", 5)
    f = lapply(strsplit(as.character(x.ligation), ''), "[", 6)

    col.matrix = cbind(a,b,c,d,e,f)
    a.nuc = sapply(1:6, function(x) sum(col.matrix[,x] == "A"))
    t.nuc = sapply(1:6, function(x) sum(col.matrix[,x] == "T"))
    c.nuc = sapply(1:6, function(x) sum(col.matrix[,x] == "C"))
    g.nuc = sapply(1:6, function(x) sum(col.matrix[,x] == "G"))
    
    cbind(a.nuc, c.nuc, g.nuc, t.nuc)
    pswm = cbind(a.nuc, c.nuc, g.nuc, t.nuc)
    outfile = file(out)
    writeLines(c("MEME version 4", "ALPHABET= ACGT", "strands: + -", " ", 
                 "Background letter frequencies (from uniform background):", 
                 "A 0.25000 C 0.25000 G 0.25000 T 0.25000", paste("MOTIF", out), " ",
                 "letter-probability matrix: alength= 4 w= 6"), outfile)
    pswm = pswm/rowSums(pswm)
    write.table(pswm, file = out, append = TRUE, quote=FALSE, row.names =FALSE, col.names = FALSE)
    system(paste('ceqlogo -i ', out, ' -m 1 > ', out, '.eps', sep=''))
}
