# ==================================
# define our source and object files
# ==================================

OBJ_DIR= ./obj

SOURCES= Main.cpp Handy.cpp
OBJECTS= Main.o Handy.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))

# ================
# compiler options
# ================

CXX= g++
CXXFLAGS= -Wall -O2
PROG= GetMersFastq
LIBS= -L/FolderWhereBamToolsIs/bamtools/lib -lbamtools -lz
LDFLAGS = -Wl,-rpath /FolderWhereBamToolsIs/bamtools/lib
INCLUDES = -I/FolderWhereBamToolsIs/bamtools/include

# ================
# build targets
# ================

.PHONY: all
all: $(OBJ_DIR) $(PROG)

$(BUILT_OBJECTS): $(SOURCES)
	@$(CXX) -c -o $@ $(*F).cpp $(CXXFLAGS) 
   
GetMersFastq: $(BUILT_OBJECTS)
	@$(CXX) $(CXXFLAGS) -o GetMersFastq $(BUILT_OBJECTS) 

$(OBJ_DIR):
	@mkdir -p $@

.PHONY: clean
clean:
	@rm -f $(OBJ_DIR)/* $(PROG)
