rule load_data:
     output:
        "data/all_matrix_list.rds"
#     conda:
#        "../envs/MOFA.yaml"
     shell:
        "Rscript scripts/mofa/load_data.R --met=../scNMT_NOMeWorkFlow/data/met --acc=../scNMT_NOMeWorkFlow/data/acc/ --rna=../scNMT_transcriptomeMapping/data/seurat/SeuratObject.rds --qcinfo=../scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --genemeta=../scNMT_transcriptomeMapping/data/gene_metadata.tsv --anno=data/anno --outdir=data"

rule run_model:
     input:
        "data/all_matrix_list.rds"
     output:
        "data/hdf5/model_trained.hdf5"
#     conda:
#        "../envs/MOFA.yaml"
     shell:
        "Rscript scripts/mofa/run.R --data={input} --name=trained --scale=FALSE --factor=10 --outdir=data"

rule analyze_model:
     input:
        "data/all_matrix_list.rds",
        "data/hdf5/model_trained.hdf5"
     output:
        "plots/full_variance.png",
	"plots/factor_variance.png",
	"plots/topweights_F1.png",
	"plots/topweights_F2.png",
	"plots/weights_HM_F1.png",
	"plots/weights_HM_F2.png",
	"plots/MOFA_UMAP.png"
#     conda:
#        "../envs/MOFA.yaml"
     shell:
        "Rscript scripts/mofa/analyze_model.R --data={input[0]} --model={input[1]} --clusters=3 --plotdir=plots"
