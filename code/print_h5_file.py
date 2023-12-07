import sys
import h5py

# ==========  PRINT .h5 FILE  ==========
# =====================================
# prints the first file given as argument
# todo: implement reading nested dataset/groups 

def main():
	if(len(sys.argv)<2):
		print("Error: No file for reading was given as parameter, exiting")
		return
	
	read_file(sys.argv[1])



def read_file(file_destination):
	print(file_destination + ":")
	h5_file = h5py.File(sys.argv[1], 'r')
	for key in h5_file.keys():
		dataset = h5_file.get(key)
		print("\t" + key + ": " + str(dataset.shape) + " of " + str(dataset.dtype) ) 
	
	h5_file.close()



if __name__ == '__main__':
	main()
