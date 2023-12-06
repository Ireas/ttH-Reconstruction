compile_my_code: delta_r_matching.cpp root2h5.py
	g++ delta_r_matching.cpp -o executable -fPIC `root-config --cflags --libs`
	./executable
	python3 root2h5.py
