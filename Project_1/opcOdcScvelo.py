import scanpy as sc
import scvelo as scv
import numpy as np
import os
from datetime import datetime

curDate=datetime.today().strftime('%Y-%m-%d')

os.chdir("/home/flyhunter/Wang/output")

combLoom=scv.read('controlStress.loom', cache=True)

adata = sc.read('opc_odc.h5ad')

adata = scv.utils.merge(adata, combLoom)


scv.pp.filter_and_normalize(adata)

scv.pp.moments(adata)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)

q1=adata.obs

scv.set_figure_params()

#adata.obs

outStream="opc_odc_allInteg_scvelo_%s.png" % curDate

outStreamCond="opc_odc_allInteg_scvelo_StrContr_%s.png" % curDate


#scv.pl.velocity_embedding_grid(adata, basis = 'umap', color= 'Annotations', figsize=(10,10), dpi=300, save='opc_odc_embedGrid.png')

scv.pl.velocity_embedding_stream(adata, basis = 'umap', color= 'Annotations', figsize=(10,10), dpi=300, save= outStream)

scv.pl.velocity_graph(adata, basis = 'umap', color= 'group', figsize=(10,10), dpi=300, save= outStreamCond)

scv.tl.recover_dynamics(adata, n_jobs = 14)
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata)

outStreamDyn="opc_odc_allInteg_scvelo_dynam_%s.png" % curDate
scv.pl.velocity_embedding_stream(adata, basis = 'umap', color= 'Annotations', figsize=(10,10), dpi=300, save= outStreamDyn)

scv.tl.latent_time(adata)
scv.pl.scatter(adata, color = 'latent_time', color_map = 'gnuplot', size = 80)

# redo PCA and umap

sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=15)


sc.tl.leiden(adata)
sc.tl.paga(adata) # cahnge to grouping algoorithm to recluster
sc.pl.paga(adata)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga')
sc.tl.leiden(adata)

sc.external.pp.bbknn(adata, batch_key='group')

sc.tl.umap(adata)
sc.pl.umap(adata, color=['group'])

# velocity simple
scv.pp.moments(adata)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.pl.velocity_embedding_stream(adata, basis = 'umap', color= 'Annotations')

# velocity scvelo
scv.tl.recover_dynamics(adata, n_jobs = 14)
scv.tl.velocity(adata, mode = 'dynamical')
scv.tl.velocity_graph(adata)


outStream="opc_odc_allInteg_scvelo_reclust_batchadj_%s.png" % curDate

outGr="opc_odc_allInteg_scvelo_reclust_StrContr_batchadj_%s.png" % curDate


scv.pl.velocity_embedding_stream(adata, basis = 'umap', color= 'Annotations', figsize=(10,10), dpi=300, save= outStream)
scv.pl.velocity_embedding_grid(adata, basis = 'umap', color= 'group', figsize=(10,10), dpi=300, save= outGr)



# clear workspace, does not work too well
for name in dir():
    if not name.startswith('_'):
        del globals()[name]

        
