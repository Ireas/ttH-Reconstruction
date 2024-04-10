match: match_jets_to_tth.cpp
	g++ match_jets_to_tth.cpp -o executable -fPIC `root-config --glibs --cflags --libs`
	./executable

convert: convert_root_to_h5.py
	python3 convert_root_to_h5.py ~/git_repos/master/output/rdataframes_output.root ~/git_repos/master/output/output.h5

plot: plot.py
	python3 plot.py ~/git_repos/master/output/rdataframes_output.root

match_convert: match_jets_to_tth.cpp convert_root_to_h5.py
	g++ match_jets_to_tth.cpp -o executable -fPIC `root-config --glibs --cflags --libs`
	./executable
	python3 convert_root_to_h5.py ~/git_repos/master/output/rdataframes_output.root ~/git_repos/master/output/output.h5

match_plot: match_jets_to_tth.cpp plot.py
	g++ match_jets_to_tth.cpp -o executable -fPIC `root-config --glibs --cflags --libs`
	./executable
	python3 plot.py ~/git_repos/master/output/rdataframes_output.root

