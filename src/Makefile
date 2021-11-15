all: spatialordinal_newmarginals_shared.so pgaspatial_newmarginals_shared.so nonspatial_newmarginals_shared.so
	# Comment here preserves the prior tab
spatialordinal_newmarginals_shared.so: spatialordinal_newmarginals_shared.cpp
	Rscript --vanilla -e "TMB::compile('spatialordinal_newmarginals_shared.cpp','-O0 -g')"
pgaspatial_newmarginals_shared.so: pgaspatial_newmarginals_shared.cpp
	Rscript --vanilla -e "TMB::compile('pgaspatial_newmarginals_shared.cpp','-O0 -g')"
nonspatial_newmarginals_shared.so: nonspatial_newmarginals_shared.cpp
	Rscript --vanilla -e "TMB::compile('nonspatial_newmarginals_shared.cpp','-O0 -g')"
clean:
	rm -rf *o
