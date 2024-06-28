match: match_jets_to_tth.cpp
	g++ match_jets_to_tth.cpp -o execute_match -fPIC `root-config --glibs --cflags --libs`
	./execute_match
	rm ./execute_match


convert: convert_root_to_h5.py
	python3 convert_root_to_h5.py ~/git_repos/master/samples/matched/_matched.root ~/git_repos/master/samples/converted/_converted.h5


inject: inject_spanet_predictions.py
	cp ~/git_repos/master/samples/matched/_matched.root ~/git_repos/master/samples/injected/_injection.root
	python3 inject_spanet_predictions.py ~/git_repos/master/samples/injected/_injection.root ~/git_repos/master/samples/predicted/_prediction.h5 ~/git_repos/master/samples/converted/_converted.h5
	

neutrino_weighting: neutrino_weighting.cpp
	g++ neutrino_weighting.cpp -o execute_neutrino_weighting -fPIC `root-config --glibs --cflags --libs`
	./execute_neutrino_weighting
	rm ./execute_neutrino_weighting



everything_temp: match_jets_to_tth.cpp convert_root_to_h5.py inject_spanet_predictions.py neutrino_weighting.cpp validate_neutrino_weighting.py
	g++ match_jets_to_tth.cpp -o execute_match -fPIC `root-config --glibs --cflags --libs`
	./execute_match
	rm ./execute_match
	python3 convert_root_to_h5.py ~/git_repos/master/samples/matched/_matched.root ~/git_repos/master/samples/converted/_converted.h5
	cp ~/git_repos/master/samples/converted/_converted.h5 ~/git_repos/master/samples/predicted/_prediction.h5
	cp ~/git_repos/master/samples/matched/_matched.root ~/git_repos/master/samples/injected/_injection.root
	python3 inject_spanet_predictions.py ~/git_repos/master/samples/injected/_injection.root ~/git_repos/master/samples/predicted/_prediction.h5 ~/git_repos/master/samples/converted/_converted.h5
	g++ neutrino_weighting.cpp -o execute_neutrino_weighting -fPIC `root-config --glibs --cflags --libs`
	./execute_neutrino_weighting
	rm ./execute_neutrino_weighting
