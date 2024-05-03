match: match_jets_to_tth.cpp
	g++ match_jets_to_tth.cpp -o execute_match -fPIC `root-config --glibs --cflags --libs`
	./execute_match

convert: convert_root_to_h5.py
	python3 convert_root_to_h5.py ~/git_repos/master/samples/matched/matched.root ~/git_repos/master/samples/converted/output.h5

match_convert: match_jets_to_tth.cpp convert_root_to_h5.py
	g++ match_jets_to_tth.cpp -o execute_match -fPIC `root-config --glibs --cflags --libs`
	./execute_match
	python3 convert_root_to_h5.py ~/git_repos/master/samples/matched/matched.root ~/git_repos/master/samples/converted/output.h5


neutrino_weighting: neutrino_weighting.cpp
	g++ neutrino_weighting.cpp -o execute_neutrino_weighting -fPIC `root-config --glibs --cflags --libs`
	./execute_neutrino_weighting
