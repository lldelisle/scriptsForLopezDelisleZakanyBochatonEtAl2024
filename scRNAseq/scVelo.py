import scvelo as scv
import numpy as np
import argparse
parser = argparse.ArgumentParser(description='Run scVelo on h5ad file and store output as table')
parser.add_argument('--input', help="input h5ad file", required=True)
parser.add_argument('--output', help="output h5ad file", required=True)
parser.add_argument('--outputCSV', help="output CSV file for metadata")
parser.add_argument('--outputVelocityUMAP', help="output CSV file for velocity on UMAP")
args = parser.parse_args()
adata = scv.read(args.input)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.tl.recover_dynamics(adata)
scv.tl.latent_time(adata)
adata.write(args.output, compression='gzip')
if args.outputCSV is not None:
  adata.obs.to_csv(args.outputCSV)
if args.outputVelocityUMAP is not None:
  np.savetxt(args.outputVelocityUMAP, adata.obsm['velocity_umap'], delimiter=",")
