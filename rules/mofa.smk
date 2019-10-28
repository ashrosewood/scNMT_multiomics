rule load_data:
     input:
        "../scNMT_transcriptomeMapping/data/SeuratObject.rds",
	"../scNMT_transcriptomeMapping/data/gene_hg19.cellRanger_metadata.tsv",
	"../test_git/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt",
        "../test_git/scNMT_NOMeWorkFlow/data/acc/body.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/acc/CGI_promoter.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/acc/CTCF.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/acc/Enhancer.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/acc/MCF7_ER_peaks.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/acc/MCF7_H3K27ac_peaks.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/acc/nonCGI_promoter.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/acc/Repressed.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/met/body.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/met/CGI_promoter.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/met/CTCF.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/met/Enhancer.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/met/MCF7_ER_peaks.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/met/MCF7_H3K27ac_peaks.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/met/nonCGI_promoter.tsv.gz",
        "../test_git/scNMT_NOMeWorkFlow/data/met/Repressed.tsv.gz"
     output:
        "data/all_matrix_list.rds"
     shell:
        "Rscript scripts/mofa/load_data.R --met=/home/groups/CEDAR/woodfin/projects/NMT-seq/MCF7_plate1_matched/scNMT_NOMeWorkFlow/data/met --acc=/home/groups/CEDAR/woodfin/projects/NMT-seq/MCF7_plate1_matched/scNMT_NOMeWorkFlow/data/acc/ --rna={input[0]} --qcinfo=/home/groups/CEDAR/woodfin/projects/NMT-seq/MCF7_plate1_matched/scNMT_NOMeWorkFlow/tables/sample_stats_qcPass.txt --genemeta={input[1]} --anno=/home/groups/CEDAR/woodfin/projects/NMT-seq/MCF7_plate1_matched/scNMT_NOMeWorkFlow/data/anno --outdir=data"

rule run_model:
     input:
        "data/all_matrix_list.rds"
     output:
        "data/hdf5/model_trained.hdf5"
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
	"plots/topweights_F3.png",
	"plots/topweights_F4.png",
	"plots/topweights_F5.png",
	"plots/topweights_F6.png",
	"plots/topweights_F7.png",
	"plots/topweights_F8.png",
	"plots/topweights_F9.png",
	"plots/topweights_F10.png",
	"plots/weights_HM_F1.png",
	"plots/weights_HM_F2.png",
	"plots/weights_HM_F3.png",
	"plots/weights_HM_F4.png",
	"plots/weights_HM_F5.png",
	"plots/weights_HM_F6.png",
	"plots/weights_HM_F7.png",
	"plots/weights_HM_F8.png",
	"plots/weights_HM_F9.png",
	"plots/weights_HM_F10.png",
	"plots/MOFA_UMAP.png"
     shell:
        "Rscript scripts/mofa/analyze_model.R --data={input[0]} --model={input[1]} --clusters=2 --plotdir=plots"
