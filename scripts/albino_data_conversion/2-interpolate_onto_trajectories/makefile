all: clean
	mkdir -p output
	make -C source
	./vnu ../1-generate_trajectory/output/trajectory.txt ../input_data/data_00400.cyl ../input_data/data_00400.lum ../input_data/opdep_2D_00400.cyl

clean:
	make clean -C source
