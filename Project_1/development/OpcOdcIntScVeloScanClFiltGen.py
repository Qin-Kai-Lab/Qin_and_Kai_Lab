import scanpy as sc
import scvelo as scv
import numpy as np
import os
from datetime import datetime
import pandas as pd

curDate=datetime.today().strftime('%Y-%m-%d')

os.chdir("/home/flyhunter/Wang/output")

combLoom=scv.read('controlStress.loom', cache=True)

umapRef = pd.read_csv('OpcOdcInt_Mon3Clust85PC_scVelMetadata.csv')

targCells = umapRef['scvelo_name'].to_list()

adata = combLoom[targCells].copy()


# add group
umapRef = umapRef.set_index('scvelo_name')

adata.obs[['group', 'cell_group']] = umapRef[['group', 'cell_group']]


# normalize data and adjust for batch effect
#scv.pp.filter_and_normalize(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pp.neighbors(adata)

sc.tl.leiden(adata)
sc.tl.paga(adata) # cahnge to grouping algoorithm to recluster
sc.pl.paga(adata)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')


sc.external.pp.bbknn(adata, batch_key='group')

sc.tl.umap(adata)
#sc.pl.umap(adata, color=['group'])

#scv.pp.moments(adata, n_pcs=85)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# velocity scvelo
scv.tl.recover_dynamics(adata, n_jobs = 14)
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata)

metadat=adata.obs

outStream="OpcOdcInt_scvelo_scanClust_filtGenes_defPC_%s.png" % curDate

outGr="OpcOdcInt_scvelo_scanClust_filtGenes_defPC_StrvsContr_%s.png" % curDate

scv.pl.velocity_embedding_stream(adata, basis = 'umap', color= 'leiden')
scv.pl.velocity_embedding_grid(adata, basis = 'umap', color= 'group')

scv.pl.velocity_embedding_stream(adata, basis = 'umap', color= 'leiden', figsize=(10,10), dpi=300, save= outStream)
scv.pl.velocity_embedding_grid(adata, basis = 'umap', color= 'group', figsize=(10,10), dpi=300, save= outGr)
