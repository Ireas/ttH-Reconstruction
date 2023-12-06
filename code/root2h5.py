import uproot
import h5py


output_file_destination = "../output/v1/output_data.root"
infile = uproot.open(output_file_destination)
print(type(infile))

keys = infile.keys()
print(keys)

events = infile['name']
print(type(events))

branches = infile['name'].keys()
for branch in branches:
	print(f"{branch:20s} {infile['name'][branch]}")
	var = events['name']
	print()
	print()
	print(var)
	print(var.array())
	print()


h5file = h5py.File('h5_file.h5', 'w')
h5file.create_dataset('test_dataset', data=events['name'].array())
h5file.close()


h5file = h5py.File('h5_file.h5', 'r')
print(h5file.keys())

for key in h5file.keys():
	print(key)
	print(type(key))
	n1 = h5file.get(key)
	print(n1)
	print(type(n1))

h5file.close()

