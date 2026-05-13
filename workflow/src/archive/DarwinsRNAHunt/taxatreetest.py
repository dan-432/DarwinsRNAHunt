from ete3 import NCBITaxa, Tree

ncbi = NCBITaxa()
tax_ids = [671071, 671068, 1729650]
tree = ncbi.get_topology(tax_ids)
tree.show()