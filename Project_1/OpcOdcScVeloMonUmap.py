import scanpy as sc
import scvelo as scv
import numpy as np
import os
from datetime import datetime
import pandas as pd

curDate=datetime.today().strftime('%Y-%m-%d')

os.chdir("/home/flyhunter/Wang/output")

combLoom=scv.read('controlStress.loom', cache=True)

X_umap = scv.load('OpcOdcInt_Mon3Clust85PC_UmapCoordScVel.csv', index_col=0)

umapRef = pd.read_csv('OpcOdcInt_Mon3Clust85PC_scVelMetadata.csv')

targCells = umapRef['scvelo_name'].to_list()

adata = combLoom[targCells].copy()

# check that cells are organized in the same order

q1 = adata.obs_names.to_list()

q1 == targCells

# add group
umapRef = umapRef.set_index('scvelo_name')

adata.obs[['group', 'cell_group']] = umapRef[['group', 'cell_group']]

w1 = adata.obs

# add monocle 3 umap to adata
adata.obsm['X_umap'] = X_umap.values

sc.pl.umap(adata, color = 'cell_group')

# normalize data and adjust for batch effect

#scv.pp.filter_and_normalize(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
sc.external.pp.bbknn(adata, batch_key='group')
#sc.tl.pca(adata)
#sc.external.pp.bbknn(adata, batch_key='group')
#scv.pp.moments(adata)
# calculate velocity
#sc.pl.umap(adata, color=['group'])
# velocity simple
#scv.pp.moments(adata)
#scv.tl.velocity(adata)
#scv.tl.velocity_graph(adata)
#scv.pl.velocity_embedding_stream(adata, basis = 'umap', color= 'cell_group')

# velocity scvelo
scv.tl.recover_dynamics(adata, n_jobs = 14)
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata)


outStream="OpcOdcInt_scvelo_monoc3Clust_filtGenes_defPC_adjGr_%s.png" % curDate

outGr="OpcOdcInt_scvelo_monoc3Clust_filtGenes_defPC_adjGr_StrvsContr_%s.png" % curDate

scv.pl.velocity_embedding_stream(adata, basis = 'umap', color= 'cell_group')
scv.pl.velocity_embedding_grid(adata, basis = 'umap', color= 'group')

scv.pl.velocity_embedding_stream(adata, basis = 'umap', color= 'cell_group', figsize=(10,10), dpi=300, save= outStream)
scv.pl.velocity_embedding_grid(adata, basis = 'umap', color= 'group', figsize=(10,10), dpi=300, save= outGr)

