import numpy as np

def entropy_between_2_genes(gene_a_vals, gene_b_vals):
  rho = np.corrcoef(gene_a_vals, gene_b_vals)[0, 1]
  entropy = -((1 + abs(rho)) / 2 * np.log2((1 + abs(rho)) / 2) + (1 - abs(rho)) / 2 * np.log2((1 - abs(rho)) / 2))
  return entropy

"""
Ho, Yen-Yi, et al. "Statistical methods for identifying differentially expressed gene combinations." Gene Function Analysis. Humana Press, 2007. 171-191.
"""
