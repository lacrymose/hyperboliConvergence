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
testCSOURCE = parallalg/algorithm/test-copy.cpp

# main() function files for running the tests for each section of the program
testCSCRIPT = parallalg/algorithm/test-copy.cpp

# main() function file for running all tests
testallCSCRIPT = test-full.cpp

# definition source files for program to be tested
progCSOURCE =


# Compiler and flags
CCMP = g++-8

COPT = -Og -D_GLIBCXX_DEBUG -fno-omit-frame-pointer# -fsanitize=address

CWARN = -Wall -Wextra -Wpedantic -Wno-unused-parameter -Wshadow

CSTD = -std=c++17 -fconcepts

LIBS = -lcppunit


#-------------------------------------
#  variable definitions

# all directories
DIRS = $(testINCLDEDIR) $(testSOURCEDIR) $(testSCRIPTDIR) $(testPROGRMDIR)

# include paths for both files to test and testing files
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
	$(CCMP) $(COPT) $(CSTD) $(CWARN) $(INCLDE) -o $@ -c $<


# each test executable depends on its own object file, the object file of the test, and all program source objects
$(testPROGRM) : $(testPROGRMDIR)%.out : $(testSCRIPTDIR)%.o $(testSOURCEDIR)%.o $(progSOURCEOBJ)
	$(CCMP) $(COPT) $(CSTD) $(CWARN) $(INCLDE) $^ $(LIBS) -o $@

# executable for all tests depends on its own object file, all the test source objects, and all program source objects
$(testallPROGRM) : $(testallSCRIPTOBJ) $(testSOURCEOBJ) $(progSOURCEOBJ)
	$(CCMP) $(COPT) $(CSTD) $(CWARN) $(INCLDE) $^ $(LIBS) -o $@


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
