import uproot
import h5py


# ==========  CONSTANTS  ==========
# =================================
ROOT_FILE_DESTINATION = "../output/v1/output_data.root"
H5_FILE_DESTINATION = "./h5_file.h5"



# ==========  READ ROOT FILE  ==========
# ======================================
root_file = uproot.open(ROOT_FILE_DESTINATION)
print("=====  ROOT FILE STRUCTURE  =====")
for tree in root_file.keys():
	print("TTree: " + tree)
	for branch in root_file[tree].keys():
		print("\tTBranch: " + branch)
		var = root_file[tree][branch]
		# print("\t" + str(var.array()))




# =========  WRITE h5 FILE  ==========
# ====================================
## Write h5 files
#h5_file = h5py.File(H5_FILE_DESTINATION, 'w')
#h5_file.create_dataset('test_dataset', data=test_events['name'].array())
#h5_file.close()
#
#
## Read h5 files
#h5_file = h5py.File(H5_FILE_DESTINATION, 'r')
#print(h5_file.keys())
#for key in h5_file.keys():
#	print(key + ": " + str(h5_file.get(key)))
#h5_file.close()
