rule annotate_acc:
     output:
        "data/acc/body.tsv.gz",
        "data/acc/CGI_promoter.tsv.gz",
        "data/acc/CTCF.tsv.gz",
        "data/acc/Enhancer.tsv.gz",
        "data/acc/MCF7_ER_peaks.tsv.gz",
        "data/acc/MCF7_H3K27ac_peaks.tsv.gz",
        "data/acc/nonCGI_promoter.tsv.gz",
        "data/acc/Repressed.tsv.gz"
     shell:
        "Rscript scripts/accmet/annotate_arw_acc.R --anno=../test_git/scNMT_NOMeWorkFlow/data/anno --raw=../test_git/scNMT_NOMeWorkFlow/bismarkSE/CX/coverage2cytosine_1based/filt/binarised --meta=../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --outdir=data/acc"

rule annotate_met:
     output:
        "data/met/body.tsv.gz",
        "data/met/CGI_promoter.tsv.gz",
        "data/met/CTCF.tsv.gz",
        "data/met/Enhancer.tsv.gz",
        "data/met/MCF7_ER_peaks.tsv.gz",
        "data/met/MCF7_H3K27ac_peaks.tsv.gz",
        "data/met/nonCGI_promoter.tsv.gz",
        "data/met/Repressed.tsv.gz"
     shell:
        "Rscript scripts/accmet/annotate_arw_met.R --anno=../test_git/scNMT_NOMeWorkFlow/data/anno --raw=../test_git/scNMT_NOMeWorkFlow/bismarkSE/CX/coverage2cytosine_1based/filt/binarised --meta=../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --outdir=data/met"

rule correlate_acc:
     input:
        "data/acc/body.tsv.gz",
        "data/acc/CGI_promoter.tsv.gz",
        "data/acc/CTCF.tsv.gz",
        "data/acc/Enhancer.tsv.gz",
        "data/acc/MCF7_ER_peaks.tsv.gz",
        "data/acc/MCF7_H3K27ac_peaks.tsv.gz",
        "data/acc/nonCGI_promoter.tsv.gz",
        "data/acc/Repressed.tsv.gz"
     output:
        "plots/cor_acc/acc_rna_correlations.pdf",
        "plots/cor_acc/acc_rna_correlations.tsv"
     shell:
        "Rscript scripts/accrna/correlations_acc.R --anno=../test_git/scNMT_NOMeWorkFlow/data/anno --meta=../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --SO=../scNMT_transcriptomeMapping/data/SeuratObject.rds --accdir=data/acc --genefile=../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv --plotdir=plots/cor_acc"

rule correlate_met:
     input:
        "data/met/body.tsv.gz",
        "data/met/CGI_promoter.tsv.gz",
        "data/met/CTCF.tsv.gz",
        "data/met/Enhancer.tsv.gz",
        "data/met/MCF7_ER_peaks.tsv.gz",
        "data/met/MCF7_H3K27ac_peaks.tsv.gz",
        "data/met/nonCGI_promoter.tsv.gz",
        "data/met/Repressed.tsv.gz"
     output:
        "plots/cor_met/met_rna_correlations.pdf",
        "plots/cor_met/met_rna_correlations.tsv"
     shell:
        "Rscript scripts/metrna/correlations_met.R --anno=../test_git/scNMT_NOMeWorkFlow/data/anno --meta=../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --SO=../scNMT_transcriptomeMapping/data/SeuratObject.rds --accdir=data/met --genefile=../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv --plotdir=plots/cor_met"

rule compare_cor:
     input:
        "plots/cor_met/met_rna_correlations.tsv",
        "plots/cor_acc/acc_rna_correlations.tsv"
     output:
        "plots/cor_accmetrna/accmetrna_correlations.scatter.pdf",
        "plots/cor_accmetrna/accmetrna_correlations.tsv",
        "plots/cor_accmetrna/accmetrna_correlations.boxplot.pdf"
     shell:
        "Rscript scripts/accmetrna/compare_corr.R --metrna={input[0]} --accrna={input[1]} --plotdir=plots/cor_accmetrna"

