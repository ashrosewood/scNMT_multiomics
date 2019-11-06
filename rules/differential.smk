rule get_groups:
     output:
        "data/groups.tsv"
     shell:
        "Rscript scripts/dif_metacc/groups_from_seurat.R --SO=../scNMT_transcriptomeMapping/data/seurat/SeuratObject.rds --out={output}"

rule dif_acc:
     input:
        "data/groups.tsv"
     output:
        "tables/difacc.tsv"
     shell:
        "Rscript scripts/dif_metacc/difacc.R --meta=../scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --acc=../scNMT_NOMeWorkFlow/data/acc --anno=data/anno --genemeta=../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv --groups={input} --out={output} --mincells=1"

rule dif_met:
     input:
        "data/groups.tsv"
     output:
        "tables/difmet.tsv"
     shell:
        "Rscript scripts/dif_metacc/difmet.R --meta=../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --met=../scNMT_NOMeWorkFlow/data/met --anno=data/anno --genemeta=../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv --groups={input} --out={output} --mincells=1"

rule plot_dif:
     input:
        "tables/difacc.tsv",
        "tables/difmet.tsv"
     output:
        "plots/acc_volcano.png",
        "plots/met_volcano.png",
        "plots/dif_scatter.png",
        "plots/dif_venn.svg"
     shell:
        "Rscript scripts/dif_metacc/plot_difmet_difacc.R --met={input[1]} --acc={input[0]}"
