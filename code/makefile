compile_and_run_all: truth_reco_matching.cpp convert_root_to_h5.py
	g++ truth_reco_matching.cpp -o executable -fPIC `root-config --cflags --libs`
	./executable
	python3 convert_root_to_h5.py

compile_and_run: truth_reco_matching.cpp
	g++ truth_reco_matching.cpp -o executable -fPIC `root-config --cflags --libs`
	./executable

compile: truth_reco_matching.cpp
	g++ truth_reco_matching.cpp -o executable -fPIC `root-config --cflags --libs`

compare: print_root_file.py print_h5_file.py
	echo ">>>>> .root file"
	python3 print_root_file.py ../output/v1/output_data.root
	echo ">>>>> .h5 file"
	python3 print_h5_file.py ../output/v1/output_data.h5
