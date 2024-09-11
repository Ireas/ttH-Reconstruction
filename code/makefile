match: match_jets_to_tth.cpp
	g++ match_jets_to_tth.cpp -o execute_match -fPIC `root-config --glibs --cflags --libs`
	./execute_match
	rm ./execute_match

neutrino_weighting: neutrino_weighting.cpp
	g++ neutrino_weighting.cpp -o execute_neutrino_weighting -fPIC `root-config --glibs --cflags --libs`
	./execute_neutrino_weighting
	rm ./execute_neutrino_weighting

all-in: match_jets_to_tth.cpp neutrino_weighting.cpp
	g++ match_jets_to_tth.cpp -o execute_match -fPIC `root-config --glibs --cflags --libs`
	./execute_match
	rm ./execute_match
	python3 merge_root.py
	python3 inject_spanet_predictions.py
	g++ neutrino_weighting.cpp -o execute_neutrino_weighting -fPIC `root-config --glibs --cflags --libs`
	./execute_neutrino_weighting
	rm ./execute_neutrino_weighting

prepare_data: prepare_data.cpp
	g++ prepare_data.cpp -o execute_prepare -fPIC `root-config --glibs --cflags --libs`
	./execute_prepare
	rm ./execute_prepare

neutrino_weighting_data: neutrino_weighting_data.cpp
	g++ neutrino_weighting_data.cpp -o execute_neutrino_weighting_data -fPIC `root-config --glibs --cflags --libs`
	./execute_neutrino_weighting_data
	rm ./execute_neutrino_weighting_data
