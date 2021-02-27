# Makefile for create_svm_vector.cpp
export PREDZINC:=${PWD}/..
export LD_LIBRARY_PATH:=${PREDZINC}/lib
CP = /bin/cp -f 

all:
	make -f myfunc.makefile
	make install -f myfunc.makefile
	make -f create_svm_vector.makefile
	make -f database_build.makefile
	make -f gating_gistPred.makefile
	make -f mtx2modm.makefile
	make -f postscan-frag-search.makefile
	make -f pssm2modm.makefile
	make -f pssm2Qij.makefile
	make -f report_zincpred.makefile
	make -f search_new.makefile
	make -f znpred-postscan.makefile
	make -f getfilepath.makefile
	make -f rootname.makefile
install:
	make install -f create_svm_vector.makefile
	make install -f database_build.makefile
	make install -f gating_gistPred.makefile
	make install -f mtx2modm.makefile
	make install -f postscan-frag-search.makefile
	make install -f pssm2modm.makefile
	make install -f pssm2Qij.makefile
	make install -f report_zincpred.makefile
	make install -f search_new.makefile
	make install -f znpred-postscan.makefile
	make install -f getfilepath.makefile
	make install -f rootname.makefile
clean:
	make clean -f myfunc.makefile
	make clean -f create_svm_vector.makefile
	make clean -f database_build.makefile
	make clean -f gating_gistPred.makefile
	make clean -f mtx2modm.makefile
	make clean -f myfunc.makefile
	make clean -f postscan-frag-search.makefile
	make clean -f pssm2modm.makefile
	make clean -f pssm2Qij.makefile
	make clean -f report_zincpred.makefile
	make clean -f search_new.makefile
	make clean -f znpred-postscan.makefile
	make clean -f getfilepath.makefile
	make clean -f rootname.makefile
