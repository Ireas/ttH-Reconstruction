import sys
import uproot # root in python

# ==========  PRINT .root FILES  ==========
# =========================================
# prints out the .root file


def main():
	if len(sys.argv)<2:
		print("Error: no file was given as argument, exiting")
		return

	print_file(sys.argv[1])


def print_file(file_destination):
	print(">> PRINTING .root FILE AT " + file_destination)
	root_file = uproot.open(sys.argv[1])
	for tree in root_file.keys():
		print("TTree: " + tree)
		for branch in root_file[tree].keys():
			content = root_file[tree][branch].array()
			print("\tTBranch: " + branch + ": " + str(content))


if __name__ == '__main__':
	main()
