match: match_jets_to_tth.cpp
	g++ match_jets_to_tth.cpp -o execute_match -fPIC `root-config --glibs --cflags --libs`
	./execute_match
	rm ./execute_match

neutrino_weighting: neutrino_weighting.cpp
	g++ neutrino_weighting.cpp -o execute_neutrino_weighting -fPIC `root-config --glibs --cflags --libs`
	./execute_neutrino_weighting
	rm ./execute_neutrino_weighting
