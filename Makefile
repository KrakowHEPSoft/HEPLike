#-----------------------------------------------------------------------------
# Definitions
# author: Marcin Chrzaszcz
DEBUG = false

# compiler
CXX = g++ -std=c++0x
#CXX = clang++ -std=c++11 # uncomment me to compile with clang (NOT TESTED!!)

# directories
BIN_DIR = bin
OBJ_DIR = obj
LIB_DIR = lib
SRC_DIR = src
INC_DIR = include
MAIN_DIR = main
TEST_DIR = test

# ROOT and RooFit
ROOTCONFIG   = root-config
ROOTCXXFLAGS = $(shell $(ROOTCONFIG) --cflags)
ROOTLIBS     = $(shell $(ROOTCONFIG) --libs) -lTreePlayer -lTMVA -lMinuit -lXMLIO -lMLP -lRIO -lRooFit -lRooStats -lRooFitCore


BOOST_LIB_OPT  =  -lboost_program_options

# GSL
GSL_INC_DIR = $(LCGGSL)/include
GSL_LIB_OPT = -L$(LCGGSL)/lib/ -lgsl -lgslcblas #-ln

# libraries and flags
LIBS = $(ROOTLIBS) $(BOOST_LIB_OPT) $(GSL_LIB_OPT)
INC = -I$(INC_DIR) -I$(BOOST_INC_DIR) -I$(GSL_INC_DIR)

# do something like 
#  make DEBUG=true
# to activate the debug options
ifeq ($(DEBUG),true)
	CXXFLAGS     = -O0 -Wall -ggdb -fPIC $(INC) $(ROOTCXXFLAGS)
else
	CXXFLAGS     = -O2  -Wall -std=c++0x -fPIC -Xlinker -zmuldefs  -Wno-unused-local-typedefs $(INC) $(ROOTCXXFLAGS)
endif

SOFLAGS      = -shared
SHLIB        = $(LIB_DIR)/libKstarLL.so


#-----------------------------------------------------------------------------
# Build lists

SRCS=$(filter-out $(wildcard $(SRC_DIR)/*/_*), $(wildcard $(SRC_DIR)/*/*.cpp))
OBJS_SRC=$(subst $(SRC_DIR), $(OBJ_DIR), $(subst .cpp,.o,$(SRCS)))
OBJS_SRC+=$(LIB_DIR)/eventdict.o # also make the TObject crap using cint
##OBJS_SRC+=$(SHLIB)
INCS_CINT=$(filter-out $(wildcard $(INC_DIR)/*/_*), $(wildcard $(INC_DIR)/pdfs/*.hpp)) # to be included in the cint comands

EXES=$(filter-out $(wildcard $(MAIN_DIR)/_*), $(wildcard $(MAIN_DIR)/*.cpp))
OBJS_EXE=$(subst $(MAIN_DIR), $(OBJ_DIR), $(subst .cpp,.o,$(EXES)))
BINS=$(subst $(MAIN_DIR), $(BIN_DIR), $(subst .cpp,,$(EXES)))

TESTS=$(filter-out $(wildcard $(TEST_DIR)_*), $(wildcard $(TEST_DIR)/*.cpp))
#OBJS_TST=$(subst $(TEST_DIR), $(OBJ_DIR), $(subst .cpp,.o,$(TESTS)))
#TESTB=$(subst $(TEST_DIR), $(BIN_DIR), $(subst .cpp,,$(TESTS)))


#-----------------------------------------------------------------------------
# Main targets
all : $(SRCS) $(OBJS_SRC) $(OBJS_EXE) $(BINS) $(SHLIB)

# connect classes, namespaces, libraries
$(SRC_DIR)/%.cc : $(INC_DIR)/%.hpp

# build all objects (classes, namespaces, libraries) i.e. everything in src
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

# $(INCS_CINT)
lib/eventdict.o : $(INCS_CINT)
	@echo -e "\nBuilt objects, now generating dictionary ..."
	rootcint -f $(LIB_DIR)/eventdict.cc -c $(INCS_CINT) include/pdfs/LinkDef.h
	$(CXX) $(CXXFLAGS) -c -I ./ -o lib/eventdict.o lib/eventdict.cc

###$(CXX) $(CXXFLAGS) -o $(LIB_DIR)/eventdict.o $(LIB_DIR)/eventdict.cc
#$(CXX) $(CXXFLAGS) -I$(INC_DIR) -o $(LIB_DIR)/eventdict.o $(LIB_DIR)/eventdict.cc
#
#cd $(LIB_DIR); $(CXX) $(CXXFLAGS) -Iinclude/ -o eventdict.o eventdict.cc
#cd $(DESTDIR); $(COMPILE) -I ../src/ -o eventdict.o eventdict.cc

# make all 'main' actions exectutables into objects
$(OBJ_DIR)/%.o : $(MAIN_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

# compile 'main' actions into executables
$(BIN_DIR)/% : $(OBJ_DIR)/%.o $(OBJS_SRC)
	$(CXX) $(CXXFLAGS) -static -o $@ -g $^ $(LIBS)


#-----------------------------------------------------------------------------
# testing code
test: $(TESTS) $(OBJS_TST) $(TESTB)

# make test exectutables into objects
$(OBJ_DIR)/%.o : $(TEST_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

# compile 'main' actions into executables
$(BIN_DIR)/% : $(OBJ_DIR)/%.o $(OBJS_SRC)
	$(CXX) $(CXXFLAGS) -o $@ -g $^ $(LIBS)

# link all objects into shared object
#$(SHLIB) : $(LIB_DIR)/%.o $(INC_DIR)/%.h

$(SHLIB) : $(OBJS_SRC) lib/eventdict.o $(INCS_CINT) lib/eventdict.h
	@echo -e "\nNow building shared object ..."
	$(CXX) $(SOFLAGS) $(CXXFLAGS) -o $(SHLIB) $(OBJS_SRC)


#-----------------------------------------------------------------------------
# to cleanup
clean:
	rm -f $(OBJ_DIR)/*.o $(OBJ_DIR)/*/*.o $(BIN_DIR)/[^C]* 
	rm -f $(LIB_DIR)/eventdict.*
	#rm -f $(SRC_DIR)/all_dict.cpp $(INC_DIR)/all_dict.h
	#rm -rf doc/html


#-----------------------------------------------------------------------------
# to help debugging makefile problems
info:
	@echo -e "ROOTCXXFLAGS:" $(ROOTCXXFLAGS) "\n"
	@echo -e "ROOTLIBS:" $(ROOTLIBS) "\n"
	@echo -e "BOOST_INC_DIR:" $(BOOST_INC_DIR) "\n"
	@echo -e "BOOST_LIB_OPT:" $(BOOST_LIB_OPT) "\n"
	@echo -e "SRCS:" $(SRCS) "\n"
	@echo -e "OBJS_SRC:" $(OBJS_SRC) "\n\n"
	@echo -e "EXES:" $(EXES) "\n"
	@echo -e "TESTS:" $(TESTS) "\n"
	@echo -e "OBJS_EXE:" $(OBJS_EXE) "\n"
	@echo -e "BINS:" $(BINS) "\n\n"
	#@echo -e "OBJS_TST:" $(OBJS_TST) "\n"
	#@echo -e "TESTB:" $(TESTB) "\n"
	@echo -e "INCS_CINT:" $(INCS_CINT) "\n"
	@echo -e $(HOST)


