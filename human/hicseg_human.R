install.packages("HiCseg")
library(HiCseg)
lines <- c("A549_NA_NA_", "HEK293_siRNA-CTCF_NA_", "HEK293_siRNA-Control_NA_", "HepG2_NA_NA_", "RAD21cv-HEK293_HRV-treated_NA_", "RAD21cv-HEK293_TEV-treated_NA_")
chroms <- c("15", "16", "17", "18", "19", "20", "21", "22", "Y")
repls <- c("1", "2")
for (line in lines) {
  for (repl in repls) {
    for (chrom in chroms) {
      matrix = as.matrix(read.table(file=paste("matrices/", line, repl, ".20000.chr", chrom, ".txt.gz", sep=""),sep="\t",header=F,stringsAsFactors=F))
      dim=dim(matrix)
      n=dim[1]
      result = HiCseg_linkC_R(n,round(n/3),"P",matrix,"D")
      write.table(result, "buff_hicseg.txt", sep="\t", quote=F)
      results<-read.table("buff_hicseg.txt",header=T,sep="\t",stringsAsFactors=F)
      options(scipen=999)
      res<-20000
      results_trim<-results[results$t_hat!=0,]
      TADs<-data.frame("start"=(results_trim$t_hat[-(nrow(results_trim))]*res),"end"=(results_trim$t_hat[-1]*res))
      write.table(TADs, file=paste("yielded/hicseg_", line, repl, ".20000.chr", chrom, ".txt", sep=""), quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE)
    }
  }
}
