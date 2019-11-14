rule correlate_acc:
     output:
        "plots/cor_acc/acc_rna_correlations.pdf",
        "plots/cor_acc/acc_rna_correlations.tsv"
     conda:
        "../envs/NMT_miltiomeCorr.yaml"
     shell:
        "Rscript scripts/accrna/correlations_acc.R --anno=../test_git/scNMT_NOMeWorkFlow/data/anno --meta=../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --SO=../scNMT_transcriptomeMapping/data/SeuratObject.rds --accdir=../test_git/scNMT_NOMeWorkFlow/data/acc --genefile=../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv --plotdir=plots/cor_acc"

rule correlate_met:
     output:
        "plots/cor_met/met_rna_correlations.pdf",
        "plots/cor_met/met_rna_correlations.tsv"
     conda:
        "../envs/NMT_miltiomeCorr.yaml"	
     shell:
        "Rscript scripts/metrna/correlations_met.R --anno=../test_git/scNMT_NOMeWorkFlow/data/anno --meta=../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --SO=../scNMT_transcriptomeMapping/data/SeuratObject.rds --accdir=../test_git/scNMT_NOMeWorkFlow/data/met --genefile=../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv --plotdir=plots/cor_met"

rule compare_cor:
     input:
        "plots/cor_met/met_rna_correlations.tsv",
        "plots/cor_acc/acc_rna_correlations.tsv"
     output:
        "plots/cor_accmetrna/accmetrna_correlations.scatter.pdf",
        "plots/cor_accmetrna/accmetrna_correlations.tsv",
        "plots/cor_accmetrna/accmetrna_correlations.boxplot.pdf"
     conda:
        "../envs/NMT_miltiomeCorr.yaml"
     shell:
        "Rscript scripts/accmetrna/compare_corr.R --metrna={input[0]} --accrna={input[1]} --plotdir=plots/cor_accmetrna"

