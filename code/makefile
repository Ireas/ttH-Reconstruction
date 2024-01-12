compile: match_jets_to_tth.cpp
	g++ match_jets_to_tth.cpp -o executable -fPIC `root-config --glibs --cflags --libs`

compile_run: match_jets_to_tth.cpp
	g++ match_jets_to_tth.cpp -o executable -fPIC `root-config --glibs --cflags --libs`
	./executable

compile_run_convert: match_jets_to_tth.cpp convert_root_to_h5.py
	g++ match_jets_to_tth.cpp -o executable -fPIC `root-config --glibs --cflags --libs`
	./executable
	python3 convert_root_to_h5.py ~/git_repos/master/output/rdataframes_output.root ~/git_repos/master/output/output.h5

convert: convert_root_to_h5.py
	python3 convert_root_to_h5.py ~/git_repos/master/output/rdataframes_output.root ~/git_repos/master/output/output.h5


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
