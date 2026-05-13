from ete3 import NCBITaxa, Tree
import sys

from DarwinsRNAHunt.downloadallfastafiles import download_file

input_file = sys.argv[1]
output_file = sys.argv[2]

ncbi = NCBITaxa()
tax_ids = [671071, 671068, 1729650]
tree = ncbi.get_topology(tax_ids)
tree.show()

tree.render(output_file)
