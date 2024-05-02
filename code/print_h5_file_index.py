import sys
import h5py

# ==========  PRINT .h5 FILE  ==========
# ======================================
# prints the first file given as argument in terminal


def main():
	if(len(sys.argv)<2):
		print("Error: no file was given as argument, exiting")
		return
	if(len(sys.argv)<3):
		print("Error: no index was given as argument, exiting")
		return
	
	print_file(  sys.argv[1], int(sys.argv[2]) )


def print_file(file_destination, index):
	print(">> PRINTING .h5 FILE AT " + file_destination)
	h5_file = h5py.File(sys.argv[1], 'r')
	print_group("FILE", h5_file, 0, index)	#treat file level as group	
	h5_file.close()


def print_group(key, group, spacing, index):
	print("|" + spacing*"-" + ">" + key)
	spacing+= 1

	for entry_key in group:
		entry = group[entry_key]
		if isinstance(entry, h5py.Group):
			print_group(entry_key, entry, spacing, index)
		elif isinstance(entry, h5py.Dataset):
			print_dataset(entry_key, entry, spacing, index)
		else:
			print("Error: key " + key + " yields no group or dataset")


def print_dataset(key, dataset, spacing, index):
		print("|" + spacing*"-" + " " + key + ": " + str(dataset[index]) ) # dtype=object is array 
	

if __name__ == '__main__':
	main()
