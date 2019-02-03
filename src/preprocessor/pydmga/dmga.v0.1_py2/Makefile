#this is a facade that redirects to src/Makefile

# Tell make that these are phony targets
.PHONY: all clean

default: all

PYTHON_FLAGS = 
# uncomment if want to use custom flags to python install script
# instalation of dmga2py is required to pydmga to work
#PYTHON_FLAGS = --prefix=~

# Build all of the executable files
all:
	$(MAKE) -C src/dmga
	@echo "*** dmga compiled sucessfuly!                                 ***"
	cd src/dmga2py/; python setup.py install $(PYTHON_FLAGS)
	@echo "*** dmga2py installed sucessfuly!                             ***"
	@echo "*** pydmga is under src/pydmga. You may want to copy it to    ***"
	@echo "*** where all your python packages are stored or to add it to ***"
	@echo "*** system path variable.                                     ***"

# Clean up the executable files
clean:
	$(MAKE) -C src/dmga clean
	
# Clean up the executable files
dmgaclean:
	$(MAKE) -C src/dmga dmgaclean
	
# Clean up the executable files
voroclean:
	$(MAKE) -C src/dmga voroclean	
