import sys
import h5py

# ==========  PRINT .h5 FILE  ==========
# ======================================
# prints the first file given as argument in terminal


def main():
	if(len(sys.argv)<2):
		print("Error: No file for reading was given as parameter, exiting")
		return
	
	print_file(sys.argv[1])


def print_file(file_destination):
	print("PRINTING FILE AT " + file_destination)
	h5_file = h5py.File(sys.argv[1], 'r')
	print_group("FILE", h5_file, 0)	#treat file level as group	
	h5_file.close()


def print_group(key, group, spacing):
	print(spacing*" " + "->" + key)
	spacing+= len(key)+2

	for entry_key in group:
		entry = group[entry_key]
		if isinstance(entry, h5py.Group):
			print_group(entry_key, entry, spacing)
		elif isinstance(entry, h5py.Dataset):
			print_dataset(entry_key, entry, spacing)
		else:
			print("Error: key " + key + " yields no group or dataset")


def print_dataset(key, dataset, spacing):
		print(spacing*" " + key + ": " + str(dataset.shape) + " of " + str(dataset.dtype) ) # dtype=object is array 
	

if __name__ == '__main__':
	main()
