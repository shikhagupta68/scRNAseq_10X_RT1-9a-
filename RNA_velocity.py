{\rtf1\ansi\ansicpg1252\cocoartf2636
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
{\*\listtable{\list\listtemplateid1\listhybrid{\listlevel\levelnfc0\levelnfcn0\leveljc0\leveljcn0\levelfollow0\levelstartat1\levelspace360\levelindent0{\*\levelmarker \{decimal\}}{\leveltext\leveltemplateid1\'01\'00;}{\levelnumbers\'01;}\fi-360\li720\lin720 }{\listname ;}\listid1}}
{\*\listoverridetable{\listoverride\listid1\listoverridecount0\ls1}}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural\partightenfactor0

\f0\fs24 \cf0 # RNA velocity analysis using velocyto.py\
\
#generating loompy file\
\
loompy fromfq new_t1_merged.loom t1 mouse_GRCh38_gencode.v31 metadata.tab SIGAB6_L001_R1_001.fastq.gz SIGAB6_L001_R2_001.fastq.gz\
 \
import velocyto as vcy\
\
t1vlm = vcy.VelocytoLoom("new_t1_merged.loom")   # initiating velocyto object\
\pard\tx220\tx720\pardeftab720\li720\fi-720\sa100\partightenfactor0
\ls1\ilvl0\cf0 \expnd0\expndtw0\kerning0
\
suffixPrefixMapping = \{"IN46_":"-3","IN42_":"-1","IN56_":"-4","IN44_":"-2"\}                                                                                                       \
tSNE1_t1 = [None for i in range(t1vlm.ca['CellID'].shape[0])]                                                                                                                                                                                                                    \
cluster_t1 = [None for i in range(t1vlm.ca['CellID'].shape[0])]                                                                                                                   \
tSNE2_t1 = [None for i in range(t1vlm.ca['CellID'].shape[0])]                                                                                                                     \
\
# reading tsne embedding file\
seurattSNE_t1 = pd.read_csv('/home/t1_embed_tsne.csv', index_col = 0)                                                                                     \
\
#reading cluster ids file\
seuratcluster_t1 = pd.read_csv('/home/t1_cluster_names.csv', index_col = 0)                                                                               \
\
prefixLength = 5\
\
#generating tSNE1, tSNE2 coordinates and cluster ids\
\
for i in range(t1vlm.ca['CellID'].shape[0]):                                                                                                                                      \
...     prefixedCellID = t1vlm.ca['CellID'][i]                                                                                                                                        \
...     currPrefix = prefixedCellID[0:prefixLength]                                                                                                                                   \
...     cellID = prefixedCellID[prefixLength:] + suffixPrefixMapping[currPrefix]                                                                                                      \
...     if cellID in seurattSNE_t1.index:                                                                                                                                             \
...             tSNE1_t1[i] = seurattSNE_t1.loc[cellID]['tSNE_1']                                                                                                                     \
...             tSNE2_t1[i] = seurattSNE_t1.loc[cellID]['tSNE_2']                                                                                                                     \
...             cluster_t1[i] = seuratcluster_t1.loc[cellID]['amlseu_500_3_t1.ident']    \
... \
\
#assigning tSNE1, tSNE2 coordinates and Cluster ids to velocyto object\
\
t1vlm.ca['TSNE1']= tSNE1_t1  \
t1vlm.ca['TSNE2']= tSNE2_t1  \
t1vlm.ca['ClusterName']= cluster_t1\
\
#pre-processing and filtering \
\
t1vlm.score_detection_levels(min_expr_counts=40, min_cells_express=30)                                                                                                            \
t1vlm.filter_genes(by_detection_levels=True)\
t1vlm.score_cv_vs_mean(3000, plot=True, max_expr_avg=35)                                                                                                                         \
t1vlm.filter_genes(by_cv_vs_mean=True)\
\
#normalizing data by size \
\
t1vlm._normalize_S(relative_size=t1vlm.S.sum(0),                                                                                                                                  \
             target_size=t1vlm.S.sum(0).mean())                                                                                                                                       \
t1vlm._normalize_U(relative_size=t1vlm.U.sum(0),                                                                                                                                  \
             target_size=t1vlm.U.sum(0).mean())  \
\
#Preparation for Gamma fit\
\
t1vlm.perform_PCA()                                                                                                                                                                     \
t1vlm.knn_imputation(n_pca_dims=20, balanced=False, n_jobs=16)                                                                                                                          \
\
#Gamma fit and extrapolation\
\
t1vlm.fit_gammas()                                                                                                                                                                      \
t1vlm.predict_U()\
\
#calculation of RNA velocity\
\
t1vlm.calculate_velocity()\
t1vlm.calculate_shift(assumption="constant_velocity")                                                                                                                                   \
t1vlm.extrapolate_cell_at_t(delta_t=1.)                                                                                                                                                 \
\
#projection of velocity onto embeddings\
\
from sklearn.manifold import TSNE                                                                                                                                                       \
bh_tsne = TSNE()                                                                                                                                                                        \
t1vlm.ts = bh_tsne.fit_transform(t1vlm.pcs[:, :25])                                                                                                                                     \
t1vlm.estimate_transition_prob(hidim="Sx_sz", embed="ts", transform="sqrt", psc=1, knn_random=True, sampled_fraction=0.5)\
t1vlm.calculate_embedding_shift(sigma_corr = 0.05, expression_scaling=True)                                                                                                             \
t1vlm.calculate_grid_arrows(smooth=0.8, steps=(40, 40), n_neighbors=300)                                                                                                                \
\
#plotting RNA velocity analysis\
\
import matplotlib.pyplot as plt                                                                                                                                                         \
plt.figure(None,(20,10))\
t1vlm.set_clusters(cluster_t1)       # colouring by cluster ids                                                                                                                                                   \
color_array =  [None for i in range(t1vlm.cluster_ix.shape[0])]                                                                                                                         \
colo_dic = \{0:"r",1:"#889E19",2:"c",3:"g",4:"#E969E5",5:"w"\} \
for i in range(t1vlm.cluster_ix.shape[0]):                                                                                                                                              \
...     color_array[i] = colo_dic[t1vlm.cluster_ix[i]]                                                                                                                                      \
... \
scatter_kwargs_dict = \{"c":color_array, "alpha":1\}                                                                                                                                      \
t1vlm.plot_grid_arrows(scatter_kwargs_dict = scatter_kwargs_dict, plot_random = False)      \
plt.savefig("t1_updated_no_random.pdf")    \
plt.close()\
}