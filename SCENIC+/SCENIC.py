#SCENIC+
import sys
import os
import pandas as pd
import scanpy as sc
import pycisTopic
import matplotlib.pyplot as plt

work_dir="/home/twyang2/scenic/"
#RNA-seq
#data load
adata = sc.read_h5ad("/home/twyang2/scenic/Robjects/adata.h5ad")
#Data Preprocessing
adata.obs['cluster.names'] = adata.obs['rna_clusters']
adata.obs['cluster.names'] = (adata.obs['cluster.names'] .astype('category'))
#
adata.raw = adata
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
#
adata.write_h5ad("/home/twyang2/scenic/data/adata.h5ad")
#ATAC-seq
cell_data=adata.obs
cell_data["celltype"]=cell_data["cluster.names"].astype(str)
# Add barcode column
cell_data['barcode'] = [x.split('_')[0] for x in cell_data.index.tolist()]
#
cell_data = cell_data.rename(cell_data["barcode"].astype(str) + "-" + cell_data["orig.ident"].astype(str))
#save cell_data as tsv file
cell_data.to_csv("/home/twyang2/scenic/data/cell_data.tsv", sep='\t', index=True)
#1.Getting pseudobulk profiles
fragments_dict = {"MCMV_D7":"/home/twyang2/cellranger/counts/MCMV_D7/outs/atac_fragments.tsv.gz","MCMV_D90":"/home/twyang2/cellranger/counts/MCMV_D90/outs/atac_fragments.tsv.gz",'LCMV_Arm_D7': "/home/twyang2/cellranger/counts/LCMV_Arm_D7/outs/atac_fragments.tsv.gz",'LCMV_Arm_D21': "/home/twyang2/cellranger/counts/LCMV_Arm_D21/outs/atac_fragments.tsv.gz","LCMV_C13_D7":"/home/twyang2/cellranger/counts/LCMV_C13_D7/outs/atac_fragments.tsv.gz","LCMV_C13_D21":"/home/twyang2/cellranger/counts/LCMV_C13_D21/outs/atac_fragments.tsv.gz"
  }
#Chr
import pyranges as pr
import requests
target_url='http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes'
chromsizes=pd.read_csv(target_url, sep='\t', header=None)
chromsizes.columns=['Chromosome', 'End']
chromsizes['Start']=[0]*chromsizes.shape[0]
chromsizes=chromsizes.loc[:,['Chromosome', 'Start', 'End']]
# Exceptionally in this case, to agree with CellRangerARC annotations
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].replace('_random', '') for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x].split('_')[1] if len(chromsizes['Chromosome'][x].split('_')) > 1 else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes['Chromosome'] = [chromsizes['Chromosome'][x]+'.1' if 'chr' not in chromsizes['Chromosome'][x] else chromsizes['Chromosome'][x] for x in range(len(chromsizes['Chromosome']))]
chromsizes=pr.PyRanges(chromsizes)
#the psuedobulk profiles
from pycisTopic.pseudobulk_peak_calling import export_pseudobulk
os.makedirs("/home/twyang2/scenic/Outputs/consensus_peak_calling")
os.makedirs("/home/twyang2/scenic/Outputs/consensus_peak_calling/pseudobulk_bed_files")
os.makedirs("/home/twyang2/scenic/Outputs/consensus_peak_calling/pseudobulk_bw_files")

bw_paths, bed_paths = export_pseudobulk(input_data = cell_data,variable = 'celltype', 
                 sample_id_col = 'orig.ident',
                 chromsizes = chromsizes,
                 bed_path = "/home/twyang2/scenic/Outputs/consensus_peak_calling/pseudobulk_bed_files",
                 bigwig_path = "/home/twyang2/scenic/Outputs/consensus_peak_calling/pseudobulk_bw_files",
                 path_to_fragments = fragments_dict,
                 n_cpu = 10,
                 normalize_bigwig = True,
                 split_pattern = '-')

#write tsv
with open("/home/twyang2/scenic/Outputs/consensus_peak_calling/bw_paths.tsv", "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")

with open("/home/twyang2/scenic/Outputs/consensus_peak_calling/bed_paths.tsv", "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")          
#
bw_paths = {}
with open("/home/twyang2/scenic/Outputs/consensus_peak_calling/bw_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bw_paths.update({v: p})
        
bed_paths = {}
with open("/home/twyang2/scenic/Outputs/consensus_peak_calling/bed_paths.tsv") as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})

##2.Inferring consensus peaks
from pycisTopic.pseudobulk_peak_calling import *
macs_path = "macs2"

os.makedirs("/home/twyang2/scenic/Outputs/consensus_peak_calling/MACS", exist_ok = True)
#Run peakcalling
narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = "/home/twyang2/scenic/Outputs/consensus_peak_calling/MACS",
    genome_size = 'mm',
    n_cpu = 19,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05
)  
#
# Other param
peak_half_width=250
path_to_blacklist="/home/twyang2/scenic/RawData/mm10-blacklist.v2.bed"
# Get consensus peaks
from pycisTopic.iterative_peak_calling import get_consensus_peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes,
    path_to_blacklist = path_to_blacklist) 
#Write to bed
consensus_peaks.to_bed(
    path = "/home/twyang2/scenic/Outputs/consensus_peak_calling/consensus_regions.bed",
    keep =True,
    compression = 'infer',
    chain = False)  
    
#3. QC
import pybiomart as pbm
#Bash "get_tss.sh"
#Run cistopic
regions_bed_filename =  "/home/twyang2/scenic/Outputs/consensus_peak_calling/consensus_regions.bed"
tss_bed_filename = "/home/twyang2/scenic/Outputs/qc/tss.bed"
pycistopic_qc_commands_filename = "pycistopic_qc_commands.txt"
# Create text file with all pycistopic qc command lines.
with open(pycistopic_qc_commands_filename, "w") as fh:
    for sample, fragment_filename in fragments_dict.items():
        print(
          "pycistopic qc",
          f"--fragments {fragment_filename}",
          f"--regions {regions_bed_filename}",
          f"--tss {tss_bed_filename}",   
          f"--output {os.path.join(work_dir, 'Outputs/qc')}/{sample}",
          sep=" ",file=fh,
        )
#Run bash "pycistopic_qc_commands.sh"   
  
#Visulation
#Sample Level
from pycisTopic.plotting.qc_plot import plot_sample_stats, plot_barcode_stats
for sample_id in fragments_dict:
    fig = plot_sample_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = os.path.join(work_dir, "Outputs/qc")
    )
    output_dir = os.path.join(work_dir, "figs/", sample_id)
    os.makedirs(output_dir, exist_ok=True)
    plt.show()
    plt.savefig(os.path.join(output_dir, "QC.png"))
    
#Barcode Level
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id = sample_id,
            pycistopic_qc_output_dir = os.path.join(work_dir, "Outputs/qc"),
            unique_fragments_threshold = None, # use automatic thresholding
            tss_enrichment_threshold = None, # use automatic thresholding
            frip_threshold = 0,
            use_automatic_thresholds = True,
    )

for sample_id in fragments_dict:
    fig = plot_barcode_stats(
        sample_id = sample_id,
        pycistopic_qc_output_dir = os.path.join(work_dir, "Outputs/qc"),
        bc_passing_filters = sample_id_to_barcodes_passing_filters[sample_id],
        detailed_title = False,
        **sample_id_to_thresholds[sample_id]
    )
    plt.show()
    plt.savefig(os.path.join(work_dir, "figs/", sample_id, "QC2.png"))
    
#4. Create cisTopic object
path_to_regions = os.path.join(work_dir, "Outputs/consensus_peak_calling/consensus_regions.bed")
path_to_blacklist = os.path.join(work_dir, "RawData/mm10-blacklist.v2.bed")
pycistopic_qc_output_dir = os.path.join(work_dir, "Outputs/qc")
#
from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import pyarrow.parquet as pq
#
cistopic_obj_list = []
for sample_id in fragments_dict:
    sample_metrics = pq.read_table(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')).to_pandas().set_index("CB").loc[ sample_id_to_barcodes_passing_filters[sample_id] ]
    valid_barcodes = sample_id_to_barcodes_passing_filters[sample_id]
    path_to_fragments =  fragments_dict[sample_id]
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = path_to_fragments,
        path_to_regions = path_to_regions,
        path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = valid_barcodes,
        n_cpu = 20,
        project = sample_id,
        split_pattern = "-"
    )
    cistopic_obj_list.append(cistopic_obj)
    
# need to merge the cistopic objects from the list
cistopic_obj = cistopic_obj_list[0]
cistopic_obj.merge(cistopic_obj_list[1:])

# check that it worked
cistopic_obj.cell_data
cistopic_obj.cell_names[0:5] 

#
import pandas as pd
cell_data = pd.read_table(os.path.join(work_dir, "data/cell_data.tsv"), index_col = 0)
cistopic_obj.add_cell_data(cell_data, split_pattern='___')
cistopic_obj.cell_data.head()
#
import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join(work_dir, "Outputs/cistopic_obj.pkl"), "wb")
)
# 3.Run cisTopic models
#
work_dir="/home/twyang2/scenic/"
import pickle
#
models = pickle.load(open(os.path.join(work_dir, "Outputs/models.pkl"), "rb"))
cistopic_obj = pickle.load(open(os.path.join(work_dir, "Outputs/cistopic_obj.pkl"), "rb"))
#
from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models,
    select_model = 200,
    return_model = True,
)
#
cistopic_obj.add_LDA_model(model)
pickle.dump(
    cistopic_obj,
    open(os.path.join(work_dir, "Outputs/cistopic_obj.pkl"), "wb")
)
#
############Clustering and Visualisation
import sys
import os
import pandas as pd
import scanpy as sc
import pycisTopic
import matplotlib.pyplot as plt
import pickle
work_dir="/home/twyang2/scenic/"
cistopic_obj = pickle.load(open(os.path.join(work_dir, "Outputs/cistopic_obj.pkl"), "rb"))
#
from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)
#
cistopic_obj.projections["cell"] = {}
#
#
run_umap(
    cistopic_obj,
    target  = 'cell',
    scale=True)
#
cistopic_obj.cell_data["celltype"]=cistopic_obj.cell_data["celltype"].astype(str)
#
plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=['celltype', 'orig.ident'],
    target='cell', num_columns=4,
    text_size=3,
    dot_size=1)
#
###Topic binarization and QC
from pycisTopic.topic_binarization import binarize_topics
#
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, num_columns=5
)
#
region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)
#
binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)
#
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img
#
topic_qc_metrics = compute_topic_metrics(cistopic_obj)
#
fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)
#
# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1
plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.savefig("figs/topic_stats.png")
plt.show()
#

topic_annot = topic_annotation(
    cistopic_obj,
    annot_var='celltype',
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)

topic_annot

## Calculating DARs
import numpy as np
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
#
variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True, save=os.path.join(work_dir,"figs/variable_features.png")
)
len(variable_regions)
#
markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable='celltype',
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir='/home/twyang2/outputs',
    split_pattern = '-'
)
#
pickle.dump(
    cistopic_obj,
    open(os.path.join(work_dir, "Outputs/cistopic_obj.pkl"), "wb")
)
#
print("Number of DARs found:")
print("---------------------")
for x in markers_dict:
    print(f"  {x}: {len(markers_dict[x])}")
    
### save region sets
os.makedirs(os.path.join(work_dir, "Outputs/region_sets"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "Outputs/region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "Outputs/region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "Outputs/region_sets", "DARs_cell_type"), exist_ok = True)
#
from pycisTopic.utils import region_names_to_coordinates
#
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(work_dir, "Outputs/region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
#
for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(work_dir, "Outputs/region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )
#
for cell_type in markers_dict:
    region_names_to_coordinates(
        markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(work_dir, "Outputs/region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )
###Create cistarget database
#create fasta(bash)
#create database(bash)
###Scenicplus pipeline
#
scenicplus
mkdir -p scplus_pipeline
scenicplus init_snakemake --out_dir scplus_pipeline
tree scplus_pipeline/
mkdir -p outs
mkdir -p tmp
#create config.taml with path
printf "%s" "$(<scplus_pipeline/Snakemake/config/config.yaml)"
cd scplus_pipeline/Snakemake
ls

#
import sys
import os
import pandas as pd
import scanpy as sc
import pycisTopic
import matplotlib.pyplot as plt
import pickle
work_dir="/home/twyang2/scenic/"
#
adata = sc.read_h5ad(os.path.join(work_dir, "data/adata.h5ad"))
cistopic_obj = pickle.load(open(os.path.join(work_dir, "Outputs/cistopic_obj.pkl"), "rb"))

adata.obs.index

cistopic_obj.projections

adata.obs.index = adata.obs_names.str.replace(r'_.*', '-', regex=True) + adata.obs["orig.ident"].astype(str)+"___" + adata.obs["orig.ident"].astype(str)

adata.obs_names

adata.write_h5ad(os.path.join(work_dir, "data/adata1.h5ad"))
#

###Run pipeline
#bash
cd scplus_pipeline/Snakemake
#Run
snakemake --cores 20

#python
import mudata
scplus_mdata = mudata.read("outs/scplusmdata.h5mu")
scplus_mdata.uns["direct_e_regulon_metadata"]
scplus_mdata.uns["extended_e_regulon_metadata"]
import scanpy as sc
import anndata
eRegulon_gene_AUC = anndata.concat(
    [scplus_mdata["direct_gene_based_AUC"], scplus_mdata["extended_gene_based_AUC"]],
    axis = 1,
)
eRegulon_gene_AUC.obs = scplus_mdata.obs.loc[eRegulon_gene_AUC.obs_names]
sc.pp.neighbors(eRegulon_gene_AUC, use_rep = "X")
sc.tl.umap(eRegulon_gene_AUC)
sc.pl.umap(eRegulon_gene_AUC, color = "scRNA_counts:celltype", save="-celltype.png")
#
# eRegulon specificity score
from scenicplus.RSS import (regulon_specificity_scores, plot_rss)
rss = regulon_specificity_scores(
    scplus_mudata = scplus_mdata,
    variable = "scRNA_counts:celltype",
    modalities = ["direct_gene_based_AUC", "extended_gene_based_AUC"]
)
plot_rss(
    data_matrix = rss,
    top_n = 3,
    num_columns = 2, fontsize = 8, save="figs/plot_rss.png"
)# Plot eRegulon enrichment scores
sc.pl.umap(eRegulon_gene_AUC, color = list(set([x for xs in [rss.loc[ct].sort_values()[0:2].index for ct in rss.index] for x in xs ])),save="-eRegulons.png")