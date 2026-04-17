rep1=read.delim('data/SummaryMisoBayesianFactorRep1.tsv')
rep2=read.delim('data/SummaryMisoBayesianFactorRep2.tsv')

r1_top5=rep1[order(rep1$MostExtremeBF, decreasing = TRUE), ][1:5, ]
r1_top5genes=sub("^gene:", "", r1_top5$event_name)
print(r1_top5genes)

r2_top5=rep2[order(rep2$MostExtremeBF, decreasing = TRUE), ][1:5, ]
r2_top5genes=sub("^gene:", "", r2_top5$event_name)
print(r2_top5genes)

for(val in r1_top5genes){
  if(val %in% r2_top5genes){
    print(val)
  }
}
