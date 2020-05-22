#-------------------------------------
# USER INPUTS:

# directories
INCLDEDIR = inc/#			header files
SOURCEDIR = src/#			class / function source files
SCRIPTDIR = script/#		main() function source files
PROGRMDIR = progrm/#		executables

# Class / function definition source files
CSOURCE = 


# main() function files
CSCRIPT = test.cpp \
	testVariableSet.cpp \
	testUtils.cpp \
	testSmooth.cpp \
	sodsExplicitEuler1.cpp \
	sodsDualTime.cpp \
	sods.cpp \
	testGenericVariableSet.cpp \
	riemannProblem.cpp \
	testScalarAdvection.cpp

# compiler
CCMP = g++-8

# debug/optimisation, warning and language standard flags
#COPT = -Og -D_GLIBCXX_DEBUG
COPT = -O3

CWARN = -Wall -Wextra -Wpedantic -Wno-unused-parameter

CSTD = -std=c++17 -fconcepts

# external libraries, eg lapack, blas, openmp ...
LIBS =


#-------------------------------------
#  variable definitions

DIRS = $(INCLDEDIR) $(SOURCEDIR) $(SCRIPTDIR) $(PROGRMDIR)

INCLDE = $(addprefix -I,$(INCLDEDIR))

# full paths for source, script and executable files
SOURCE = $(addprefix $(SOURCEDIR),$(CSOURCE))
SCRIPT = $(addprefix $(SCRIPTDIR),$(CSCRIPT))
PROGRM = $(addprefix $(PROGRMDIR),$(CSCRIPT:.cpp=.out))

# names of scripts (no suffix)
PNAMES = $(CSCRIPT:.cpp=)

# object files
SOURCEOBJ = $(SOURCE:.cpp=.o)
SCRIPTOBJ = $(SCRIPT:.cpp=.o)

OBJS = $(SOURCEOBJ) $(SCRIPTOBJ)


#-------------------------------------
# compilation recipes

# default
%.o : %.cpp
	$(CCMP) $(COPT) $(CSTD) $(CWARN) $(INCLDE) -o $@ -c $<

# each executable depends on its own object file, and all source objects

$(PROGRM) : $(PROGRMDIR)%.out : $(SCRIPTDIR)%.o $(SOURCEOBJ)
	$(CCMP) $(COPT) $(CSTD) $(CWARN) $(INCLDE) -o $@ $^ $(LIBS)


#-------------------------------------
# misc recipes

.PHONY : all clean echo obj mkdir names $(PNAMES)

# make <pname> will compile only the executable 'pname.out'
$(PNAMES) : % : $(PROGRMDIR)%.out

# make all executables
all: $(PROGRM)

# make source objects
obj: $(SOURCEOBJ)

# create required directories
mkdir:
	mkdir $(DIRS)

# print names of all executables to standard output
names:
	@for name in $(PNAMES); do echo $$name; done

# delete all non-source files
clean:
	rm -f $(OBJS) $(PROGRM)

# clean-all: clean test-clean

#-------------------------------------
# testing recipes

# .PHONY : test-%
# 
# args= `arg="$(filter-out $@,$(MAKECMDGOALS))" && echo $${arg:-${1}}`
# 
# %:
# 	@:
# 
# test-clean:
# 	$(MAKE) -C tests/ clean
# 
# test-%:
# 	$(MAKE) -C tests/ $@
# 
# test-run-%:
# 	$(MAKE) -C tests/ $(@:test-run-%=test-%)
# 	./tests/progrm/$(@:test-run-%=test-%).out
# 
# test:
# 	@for t in $(call args); do $(MAKE) -C tests/ test-$$t; done
# 
# test-run:
# 	@for t in $(call args); \
# 		do $(MAKE) -C tests/ test-$$t; \
# 		./tests/progrm/test-$$t.out; \
# 	done
# 
# echo:
# 	@echo $(call args)
# 
# #-------------------------------------
