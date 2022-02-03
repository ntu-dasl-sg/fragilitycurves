all: spatialordinal_newmarginals_shared.so pgaspatial_newmarginals_shared.so nonspatial_newmarginals_shared.so
	# Comment here preserves the prior tab
spatialordinal_newmarginals_shared.so: spatialordinal_newmarginals_shared.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('spatialordinal_newmarginals_shared.cpp','-O0 -g',safebounds=FALSE,safeunload=FALSE)"
pgaspatial_newmarginals_shared.so: pgaspatial_newmarginals_shared.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('pgaspatial_newmarginals_shared.cpp','-O0 -g',safebounds=FALSE,safeunload=FALSE)"
nonspatial_newmarginals_shared.so: nonspatial_newmarginals_shared.cpp
	$(R_HOME)/bin/Rscript --vanilla -e "TMB::compile('nonspatial_newmarginals_shared.cpp','-O0 -g',safebounds=FALSE,safeunload=FALSE)"
clean:
	rm -rf *o
