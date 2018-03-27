all:
	make -C 1-generate_trajectory
	make -C 2-interpolate_onto_trajectories
	make -C 3-separate_files_for_SQA

clean:
	make clean -C 1-generate_trajectory
	make clean -C 2-interpolate_onto_trajectories
	make clean -C 3-separate_files_for_SQA
