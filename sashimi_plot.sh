module load misopy
salloc -N1 -t2:00:00 --ntasks=8

____________________________________________________________
#MISO (comparison mode) for normal/tumor_rep2
compare_miso --compare-samples $TUMOR2 $NORMAL2 $OUTDIR2
cd $OUTDIR2/HCC1395_tumorPOS_rep2_vs_HCC1395_normalPOS_rep2/bayes-factors

##Downstream Analysis
module load r
R

miso_df = read.csv("HCC1395_tumorPOS_rep2_vs_HCC1395_normalPOS_rep2.miso_bf", header=T, sep="\t")
head(miso_df, 1)

miso_df$MostExtremeBF = sapply(as.character(miso_df$bayes_factor), function(x){
  bayes_factors = as.numeric(unlist(strsplit(x, ",", T)))
  return(max(bayes_factors))
})

miso_df$MostExtremeDiff = sapply(as.character(miso_df$diff), function(x){
  psis = as.numeric(unlist(strsplit(x, ",", T)))
  return(max(abs(psis)))
})

miso_df = miso_df[order(miso_df$MostExtremeBF, decreasing=TRUE), ]
miso_df_select = miso_df[, c("event_name", "diff", "MostExtremeDiff", "MostExtremeBF")]
write.table(miso_df_select, file="SummaryMisoBayesianFactorRep2.tsv", quote=F, row.names=F, sep="\t")

quit()

____________________________________________________________
#Rep1 Sashimi Plot
cd $HOME/isoformAnalysis
nano miso_setting_rep1.txt
'[data]
bam_prefix=/BAMS_FOR_MISO
miso_prefix=/isoformAnalysis/misoRuns
bam_files = ["HCC1395_tumorPOS_rep1.bam", "HCC1395_normalPOS_rep1.bam"]
miso_files = ["HCC1395_tumorPOS_rep1", "HCC1395_normalPOS_rep1"]
[plotting]
fig_width = 10
fig_height = 8
intron_scale = 50
exon_scale = 4
logged = True
font_size = 6.5
bar_posteriors = False
ymax = 150
nyticks = 3
nxticks = 4
show_ylabel = True
show_xlabel = True
show_posteriors = True
resolution = 0.5
posterior_bins = 40
gene_posterior_ratio = 5
colors = ["#FE9A2E", "#08088A"]
coverages = [100000, 110000]
bar_color = "b"
bf_thresholds = [0, 1, 2, 5, 10, 20]'

MISOINDEX=$HOME/isoformAnalysis/datasets/misoIndex
SETTINGFILE1=$HOME/isoformAnalysis/miso_setting_rep1.txt
OUTDIR=$HOME/isoformAnalysis/sashimiPlots

sashimi_plot --plot-event "gene:ENSG00000185721" $MISOINDEX $SETTINGFILE1 --output-dir $HOME/isoformAnalysis/sashimiPlots

cd $HOME/isoformAnalysis/sashimiPlots
mv "gene:ENSG00000185721.pdf" "rep1_ENSG00000185721.pdf"

____________________________________________________________
#Rep2 Sashimi Plot
cd $HOME/isoformAnalysis

nano miso_setting_rep2.txt
'[data]
bam_prefix=/BAMS_FOR_MISO
miso_prefix=/isoformAnalysis/misoRuns
bam_files = ["HCC1395_tumorPOS_rep2.bam", "HCC1395_normalPOS_rep2.bam"]
miso_files = ["HCC1395_tumorPOS_rep2", "HCC1395_normalPOS_rep2"]
[plotting]
fig_width = 10
fig_height = 8
intron_scale = 50
exon_scale = 4
logged = True
font_size = 6.5
bar_posteriors = False
ymax = 150
nyticks = 3
nxticks = 4
show_ylabel = True
show_xlabel = True
show_posteriors = True
resolution = 0.5
posterior_bins = 40
gene_posterior_ratio = 5
colors = ["#FE9A2E", "#08088A"]
coverages = [100000, 110000]
bar_color = "b"
bf_thresholds = [0, 1, 2, 5, 10, 20]'

MISOINDEX=$HOME/isoformAnalysis/datasets/misoIndex
SETTINGFILE2=$HOME/isoformAnalysis/miso_setting_rep2.txt
OUTDIR=$HOME/isoformAnalysis/sashimiPlots

sashimi_plot --plot-event "gene:ENSG00000185721" $MISOINDEX $SETTINGFILE2 --output-dir $OUTDIR

cd $HOME/isoformAnalysis/sashimiPlots
mv "gene:ENSG00000185721.pdf" "rep2_ENSG00000185721.pdf"

____________________________________________________________
#All Sashimi Plots (R1 + R2)
cd $HOME/isoformAnalysis
nano miso_setting_combine.txt
'[data]
bam_prefix=/BAMS_FOR_MISO
miso_prefix=/isoformAnalysis/misoRuns
bam_files = ["HCC1395_tumorPOS_rep1.bam", "HCC1395_normalPOS_rep1.bam", "HCC1395_tumorPOS_rep2.bam", "HCC1395_normalPOS_rep2.bam"]
miso_files = ["HCC1395_tumorPOS_rep1", "HCC1395_normalPOS_rep1", "HCC1395_tumorPOS_rep2", "HCC1395_normalPOS_rep2"]
[plotting]
fig_width = 10
fig_height = 8
intron_scale = 50
exon_scale = 4
logged = True
font_size = 6.5
bar_posteriors = False
ymax = 150
nyticks = 3
nxticks = 4
show_ylabel = True
show_xlabel = True
show_posteriors = True
resolution = 0.5
posterior_bins = 40
gene_posterior_ratio = 5
colors = ["#FE9A2E", "#08088A", "#FE9A2E", "#08088A"]
coverages = [100000, 110000, 100000, 110000]
bar_color = "b"
bf_thresholds = [0, 1, 2, 5, 10, 20]'

MISOINDEX=$HOME/isoformAnalysis/datasets/misoIndex
SETTINGFILEC=$HOME/isoformAnalysis/miso_setting_combine.txt
OUTDIR=$HOME/isoformAnalysis/sashimiPlots

sashimi_plot --plot-event "gene:ENSG00000185721" $MISOINDEX $SETTINGFILEC --output-dir $HOME/isoformAnalysis/sashimiPlots

cd $HOME/isoformAnalysis/sashimiPlots
mv "gene:ENSG00000185721.pdf" "ENSG00000185721.pdf"

____________________________________________________________
#All Sashimi Plots Zoomed (R1 + R2)
cd $HOME/isoformAnalysis
nano miso_setting_combine_zoom.txt
'[data]
bam_prefix=/BAMS_FOR_MISO
miso_prefix=/isoformAnalysis/misoRuns
bam_files = ["HCC1395_tumorPOS_rep1.bam", "HCC1395_normalPOS_rep1.bam", "HCC1395_tumorPOS_rep2.bam", "HCC1395_normalPOS_rep2.bam"]
miso_files = ["HCC1395_tumorPOS_rep1", "HCC1395_normalPOS_rep1", "HCC1395_tumorPOS_rep2", "HCC1395_normalPOS_rep2"]
[plotting]
fig_width = 10
fig_height = 8
intron_scale = 50
exon_scale = 4
logged = True
font_size = 6.5
bar_posteriors = False
ymax = 100
nyticks = 3
nxticks = 4
show_ylabel = True
show_xlabel = True
show_posteriors = True
resolution = 0.5
posterior_bins = 40
gene_posterior_ratio = 5
colors = ["#FE9A2E", "#08088A", "#FE9A2E", "#08088A"]
coverages = [100000, 110000, 100000, 110000]
bar_color = "b"
bf_thresholds = [0, 1, 2, 5, 10, 20]'

MISOINDEX=$HOME/isoformAnalysis/datasets/misoIndex
SETTINGFILECZ=$HOME/isoformAnalysis/miso_setting_combine_zoom.txt
OUTDIR=$HOME/isoformAnalysis/sashimiPlots

sashimi_plot --plot-event "gene:ENSG00000185721" $MISOINDEX $SETTINGFILECZ --output-dir $HOME/isoformAnalysis/sashimiPlots

cd $HOME/isoformAnalysis/sashimiPlots
mv "gene:ENSG00000185721.pdf" "ENSG00000185721_zoom.pdf"
