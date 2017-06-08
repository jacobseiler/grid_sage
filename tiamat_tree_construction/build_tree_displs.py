import h5py
import numpy as np

from mpi4py import MPI

comm= MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def write_binary_file(filename, array):
	output_file = open(filename, 'wb')
	array.tofile(output_file)
	output_file.close()

	print "Saved to file %s" %(filename)	


if __name__ == '__main__':

    filepath = '/home/msinha/scratch/tao/data_products/output/meraxes/tiamat/Tiamat_meraxes_fixed_merge_intosnapnum.h5'

    with h5py.File(filepath, 'r') as f:

	new_counts = np.empty((f['tree_counts'].shape), dtype = np.int64)
	new_displs = np.empty((f['tree_displs'].shape), dtype = np.int64)
	new_displs[0] = 0

	print "There are %d total trees." %(f['tree_counts'].shape)

	first_tree = 0
	last_tree = 50
	#for tree in xrange(0, f['tree_counts'].shape):
	for tree in xrange(first_tree + rank, last_tree + 1, size):
		print "I am rank %d and I am on tree %d" %(rank, tree)	

		gal_idx = np.arange(f['tree_displs'][tree], f['tree_displs'][tree+1], dtype = np.int64)		
		nonghost_count = len(np.where(f['galaxies']['GhostFlag'][gal_idx] == 0)[0])

		new_counts[tree] = nonghost_count

    if (rank == 0):
	counts_final = np.zeros((len(new_counts)), dtype = np.int64)	
    else:
	counts_final = None
			
    comm.Reduce([new_counts, MPI.INT], [counts_final, MPI.INT], op = MPI.SUM, root = 0)

    if (rank == 0):
	for tree in xrange(first_tree, last_tree, size):
		new_displs[tree+1] = new_displs[tree] + counts_final[tree] 
		print counts_final
		print new_displs


	filename_counts = '/lustre/projects/p004_swin/jseiler/tiamat_tree_construction/nonghost_tree_counts'
	write_binary_file(filename_counts, counts_final)

	filename_displs = '/lustre/projects/p004_swin/jseiler/tiamat_tree_construction/nonghost_tree_displs'	
	write_binary_file(filename_displs, new_displs)
