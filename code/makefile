compile: match_jets_to_tth.cpp
	g++ match_jets_to_tth.cpp -o executable -fPIC `root-config --glibs --cflags --libs`

compile_and_run: match_jets_to_tth.cpp
	g++ match_jets_to_tth.cpp -o executable -fPIC `root-config --glibs --cflags --libs`
	./executable


#old_compile_and_run_all: truth_reco_matching.cpp convert_root_to_h5.py
#	g++ truth_reco_matching.cpp -o executable -fPIC `root-config --glibs --cflags --libs`
#	./executable
#	python3 convert_root_to_h5.py
#
#old_compile_and_run: truth_reco_matching.cpp
#	g++ truth_reco_matching.cpp -o executable -fPIC `root-config --glibs --cflags --libs`
#	./executable
#
#old_compile: truth_reco_matching.cpp
#	g++ truth_reco_matching.cpp -o executable -fPIC `root-config --glibs --cflags --libs`
#
#old_compare: print_root_file.py print_h5_file.py
#	echo "\n>>>>> .root file"
#	python3 print_root_file.py ../output/v1/output_data.root
#	echo "\n>>>>> .h5 file"
#	python3 print_h5_file.py ../output/v1/output_data.h5
