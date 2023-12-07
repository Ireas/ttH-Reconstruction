import uproot
import h5py

# ==========  CONVERT .root TO .h5 FILES  ==========
# ==================================================
# prints out the .root file at fixed destination and converts it to a flat .h5 file
# todo: implement creating nested dataset/groups 

ROOT_FILE_DESTINATION = "../output/v1/output_data.root"
H5_FILE_DESTINATION = "./h5_file.h5"



def main():
	root_file = uproot.open(ROOT_FILE_DESTINATION)
	print_root_file(root_file)
	convert_to_h5(root_file)



def print_root_file(root_file):
	print("=====  ROOT FILE STRUCTURE  =====")
	for tree in root_file.keys():
		print("TTree: " + tree)
		for branch in root_file[tree].keys():
			print("\tTBranch: " + branch)
			var = root_file[tree][branch]
			# print("\t" + str(var.array()))



def convert_to_h5(root_file):
	h5_file = h5py.File(H5_FILE_DESTINATION, 'w')
	
	for tree in root_file.keys():
		for branch in root_file[tree].keys():
			var = root_file[tree][branch]
			h5_file.create_dataset(branch, data=root_file[tree][branch].array())

	h5_file.close()



if __name__ == '__main__':
	main()
