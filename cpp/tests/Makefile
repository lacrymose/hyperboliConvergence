#-------------------------------------
# USER INPUTS:

# directories
progINCLDEDIR = ../inc/#	   header files of program to test
progSOURCEDIR = ../src/#	   class / function source files of program to test

testINCLDEDIR = inc/#			test header files
testSOURCEDIR = src/#			test class / function source files

testSCRIPTDIR = script/#		test main() function source files
testPROGRMDIR = progrm/#		test executables


# definition source files for the tests for each section of the program
testCSOURCE = idealGas2D/test-Species.cpp \
              idealGas2D/test-VariableSet.cpp \
              idealGas2D/test-VariableDelta.cpp \
              idealGas2D/test-State.cpp

# main() function files for running the tests for each section of the program
testCSCRIPT = idealGas2D/test-Species.cpp \
              idealGas2D/test-VariableSet.cpp \
              idealGas2D/test-VariableDelta.cpp \
              idealGas2D/test-State.cpp

# main() function file for running all tests
testallCSCRIPT = test-full.cpp

# definition source files for program to be tested
progCSOURCE = idealGas2D/species.cpp


# Compiler and flags
CCMP = g++

COPT = -g -Wall

LIBS = -lcppunit


#-------------------------------------
#  variable definitions

DIRS = $(testINCLDEDIR) $(testSOURCEDIR) $(testSCRIPTDIR) $(testPROGRMDIR)

INCLDE = $(addprefix -I,$(testINCLDEDIR)) $(addprefix -I,$(progINCLDEDIR))

# full paths for source, script and executable files
progSOURCE = $(addprefix $(progSOURCEDIR),$(progCSOURCE))
testSOURCE = $(addprefix $(testSOURCEDIR),$(testCSOURCE))

testSCRIPT = $(addprefix $(testSCRIPTDIR),$(testCSCRIPT))
testPROGRM = $(addprefix $(testPROGRMDIR),$(testCSCRIPT:.cpp=.out))

testallSCRIPT = $(addprefix $(testSCRIPTDIR),$(testallCSCRIPT))
testallPROGRM = $(addprefix $(testPROGRMDIR),$(testallCSCRIPT:.cpp=.out))

# names of scripts (no suffix)
testPNAMES = $(testCSCRIPT:.cpp=) $(testallCSCRIPT:.cpp=)

# object files
progSOURCEOBJ = $(progSOURCE:.cpp=.o)
testSOURCEOBJ = $(testSOURCE:.cpp=.o)

testSCRIPTOBJ = $(testSCRIPT:.cpp=.o)
testallSCRIPTOBJ = $(testallSCRIPT:.cpp=.o)

OBJS = $(progSOURCEOBJ) $(testSOURCEOBJ) $(testSCRIPTOBJ) $(testallSCRIPTOBJ)


#-------------------------------------
# compilation recipes

# default
%.o : %.cpp
	$(CCMP) $(COPT) $(INCLDE) -o $@ -c $<

## wraps all program source object files into one
#allprog.o : $(progSOURCEOBJ)
#	ld -r $^ -o $@
#
## wraps all test source object files into one
#alltest.o : $(testSOURCEOBJ)
#	ld -r $^ -o $@

# each test executable depends on its own object file, the object file of the test, and all program source objects
$(testPROGRM) : $(testPROGRMDIR)%.out : $(testSCRIPTDIR)%.o $(testSOURCEDIR)%.o $(progSOURCEOBJ)
	$(CCMP) $(COPT) $(INCLDE) $^ $(LIBS) -o $@

# executable for all tests depends on its own object file, all the test source objects, and all program source objects
$(testallPROGRM) : $(testallSCRIPTOBJ) $(testSOURCEOBJ) $(progSOURCEOBJ)
	$(CCMP) $(COPT) $(INCLDE) $^ $(LIBS) -o $@


#-------------------------------------
# misc recipes

.PHONY : all clean echo obj mkdir names $(PNAMES)

# make <pname> will compile only the executable 'pname.out'
$(testPNAMES) : % : $(testPROGRMDIR)%.out

# make all executables
all: $(testPNAMES)

# make source objects
obj: $(testSOURCEOBJ) $(testallSOURCEOBJ) $(progSOURCEOBJ)

# create required directories
mkdir:
	mkdir $(DIRS)

# print names of all executables to standard output
names:
	@for name in $(testPNAMES); do echo $$name; done

# delete all non-source files
clean:
	rm -f $(OBJS) $(testPROGRM) $(testallPROGRM)

echo:
	echo

#-------------------------------------
