rule get_groups:
     input:
        "../scNMT_transcriptomeMapping/data/SeuratObject.rds"
     output:
        "data/groups.tsv"
     conda:
        "../envs/NMT_miltiomeDiff.yaml"	
     shell:
        "Rscript scripts/dif_metacc/groups_from_seurat.R --SO={input} --out={output}"

rule dif_acc:
     input:
        "data/acc",
        "data/groups.tsv"
     output:
        "tables/difacc.tsv"
     conda:
        "../envs/NMT_miltiomeDiff.yaml"		
     shell:
        "Rscript scripts/dif_metacc/difacc.R --meta=../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --acc={input[0]} --anno=../test_git/scNMT_NOMeWorkFlow/data/anno --genemeta=../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv --groups={input[1]} --out={output} --mincells=1"

rule dif_met:
     input:
        "data/met",
        "data/groups.tsv"
     output:
        "tables/difmet.tsv"
     conda:
        "../envs/NMT_miltiomeDiff.yaml"		
     shell:
        "Rscript scripts/dif_metacc/difmet.R --meta=../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --met={input[0]} --anno=../test_git/scNMT_NOMeWorkFlow/data/anno --genemeta=../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv --groups={input[1]} --out={output} --mincells=1"

rule plot_dif:
     input:
        "tables/difacc.tsv",
        "tables/difmet.tsv"
     output:
        "plots/acc_volcano.png",
        "plots/met_volcano.png",
        "plots/dif_scatter.png",
        "plots/dif_venn.svg"
     conda:
        "../envs/NMT_miltiomeDiff.yaml"		
     shell:
        "Rscript scripts/dif_metacc/plot_difmet_difacc.R --met={input[1]} --acc={input[0]}"
