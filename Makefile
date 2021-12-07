# This file manages the compilation of the Euler code. 
# All lines start with "#" are comment lines not considered 
# by the operation system

SaveName = $(USER)SaveSrc
SaveName = SaveSrc

OUT = output.o output_hg.o output_mat.o

OBS = common_block.o euler.o apply_bconds.o  set_others.o check_conv.o set_timestep.o check_grid.o smooth.o crude_guess.o flow_guess.o read_data.o  sum_fluxes.o  generate_grid.o  set_fluxes.o $(OUT) 

# Choice of comipler
FC = gfortran

# Option of Compilation

# debug option
FFLAGS = -g

# select the optimising option by delete # in the following line
#
#FFLAGS = -O2

# compile the executable and name it "Euler" using "-o"
euler : $(OBS)
	$(FC)  $(FFLAGS) -o Euler  $(OBS)  

# remove all .o files and Euler for a new start
clean:
	rm -fr *.o  Euler Save*

# save *.f90, *flow *geom in directory SaveSrc and compress it to 
# file SaveSrc.tar.gz
save:
	if [ ! -d $(SaveName) ]; then mkdir $(SaveName); fi
	cp *.f90 $(SaveName)/
	cp -r python $(SaveName)/
	cp *flow $(SaveName)/
	cp *geom $(SaveName)/
	cp Makefile $(SaveName)/
	tar -cf $(SaveName).tar $(SaveName)
	gzip $(SaveName).tar
	rm -r -f $(SaveName)

# extract files from SaveSrc.tar.gz and save them to SaveSrc
extract:
	gunzip $(SaveName).tar.gz
#tar -xvf $(SaveName).tar
	tar -xf $(SaveName).tar
# compile .f90 to .o
%.o:%.f90
	$(FC) $(FFLAGS) -c $^ -o $@

# Set flow and geom files.
0:
	cp test0_flow flow
	cp test0_geom geom

1:
	cp test1_flow flow
	cp test1_geom geom

2:
	cp test2_flow flow
	cp test2_geom geom

3:
	cp test3_flow flow
	cp test3_geom geom

4:
	cp test4_flow flow
	cp test4_geom geom

5:
	cp test5_flow flow
	cp test5_geom geom

straight:	
	cp test_straight_flow flow
	cp test_straight_geom geom