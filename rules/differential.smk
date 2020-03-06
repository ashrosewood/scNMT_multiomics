rule get_groups:
     output:
        "data/groups.tsv"
     conda:
        "../envs/NMT_miltiomeDiff.yaml"	
     shell:
        "Rscript scripts/dif_metacc/groups_from_seurat.R --SO=../scNMT_transcriptomeMapping/data/seurat/SeuratObject.rds --out={output}"

rule dif_acc:
     input:
        "data/groups.tsv"
     output:
        "tables/difaccA_vs_B.tsv"
     conda:
        "../envs/NMT_miltiomeDiff.yaml"		
     shell:
        "Rscript scripts/dif_metacc/difacc.R --meta=../scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --acc=../scNMT_NOMeWorkFlow/data/acc --anno=data/anno --genemeta=../scNMT_transcriptomeMapping/data/gene_metadata.tsv --groups={input} --out=tables/difacc --mincells=1"

rule dif_met:
     input:
        "data/groups.tsv"
     output:
        "tables/difmetA_vs_B.tsv"
     conda:
        "../envs/NMT_miltiomeDiff.yaml"		
     shell:
        "Rscript scripts/dif_metacc/difmet.R --meta=../scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --met=../scNMT_NOMeWorkFlow/data/met --anno=data/anno --genemeta=../scNMT_transcriptomeMapping/data/gene_metadata.tsv --groups={input} --out=tables/difmet --mincells=1"

rule plot_dif:
     input:
        "tables/difaccA_vs_B.tsv",
        "tables/difmetA_vs_B.tsv"
     output:
        "plots/acc_volcano_A_vs_B.png",
        "plots/met_volcano_A_vs_B.png",
        "plots/dif_scatter_A_vs_B.png",
        "plots/dif_venn_A_vs_B.svg"
     conda:
        "../envs/NMT_miltiomeDiff.yaml"	
     shell:
        "Rscript scripts/dif_metacc/plot_difmet_difacc.R --met={input[1]} --acc={input[0]}"

rule zoom:
     input:
        "data/groups.tsv",
        "tables/difaccA_vs_B.tsv",
        "tables/difmetA_vs_B.tsv"
     output:
        "plots/zoom/zoom_plot.pdf"
     conda:
        "../envs/zoom.yaml"
     shell:
        "Rscript scripts/zoom_plot.R --met=../scNMT_NOMeWorkFlow/bismarkSE/CX/coverage2cytosine_1based/filt/binarised --acc=../scNMT_NOMeWorkFlow/bismarkSE/CX/coverage2cytosine_1based/filt/binarised --rna=../scNMT_transcriptomeMapping/data/seurat/SeuratObject.rds --qcinfo=../scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --genemeta=../scNMT_transcriptomeMapping/data/gene_metadata.tsv --anno=data/anno --groups={input} --outdir=plots/zoom" 
