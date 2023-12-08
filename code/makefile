compile_and_run_all: delta_r_matching.cpp root_to_h5_conversion.py
	g++ delta_r_matching.cpp -o executable -fPIC `root-config --cflags --libs`
	./executable
	python3 root_to_h5_conversion.py

compile_and_run: delta_r_matching.cpp
	g++ delta_r_matching.cpp -o executable -fPIC `root-config --cflags --libs`
	./executable

compile: delta_r_matching.cpp
	g++ delta_r_matching.cpp -o executable -fPIC `root-config --cflags --libs`
