import scanpy as sc
import scvelo as scv
import numpy as np
import loompy


ldataC=scv.read('/home/flyhunter/Wang/data/control/velocyto/control.loom', cache=True)

adata = sc.read('/home/flyhunter/Wang/output/control.h5ad')

adata = scv.utils.merge(adata, ldata)

scv.pp.filter_and_normalize(adata)

scv.pp.moments(adata)

scv.tl.velocity(adata)

scv.tl.velocity_graph(adata)

q1=adata.obs

scv.set_figure_params()

adata.obs




scv.pl.velocity_embedding(adata, basis = 'umap', color= 'Annotations', figsize=(10,10), dpi=300, save='control.pdf')


for name in dir():
    if not name.startswith('_'):
        del globals()[name]
        

### check stress loom object

ldataS=scv.read('/home/flyhunter/Wang/data/stress/velocyto/stress.loom', cache=True)

ldata.obs_names

ldataC.var_names

ldata = ldataC.concatenate([ldataS])

loompy.combine(['/home/flyhunter/Wang/data/control/velocyto/control.loom', '/home/flyhunter/Wang/data/stress/velocyto/stress.loom'], key="Accession", output_file='controlStress.loom')


# try with combined loom file
combLoom=scv.read('controlStress.loom', cache=True)

adata = scv.utils.merge(adata, combLoom)

