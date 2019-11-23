include ${PETSC_DIR}/lib/petsc/conf/variables
#include ${PETSC_DIR}/lib/petsc/conf/rules
LIBS 	 :=  ../CSparse/Lib/libcsparse.a -lrt # -lg2c -mkl
CXXFLAGS := -O0 -g -Wall -fstrict-aliasing -fshort-enums -finline-limit=10000000000 --param large-function-growth=10000000000 --param max-inline-insns-single=10000000000 -DMPICH_IGNORE_CXX_SEEK
CXXFLAGSDEBUG := -O0 -g -Wall -Wstrict-aliasing -fstrict-aliasing -fshort-enums -finline-limit=10000000000 --param large-function-growth=10000000000 --param max-inline-insns-single=10000000000 -DMPICH_IGNORE_CXX_SEEK

CXX      := mpicxx
EXT      := cpp
TARGET   := ../bin/mspai-1.2
SOURCES  := $(wildcard *.$(EXT))
OBJECTS  := $(patsubst %.$(EXT), %.o, $(SOURCES))
IMP      := $(wildcard *.imp)
DEP			 := $(OBJECTS:.o=.d)	#one dependency file for each source



.PHONY: all
all: $(TARGET)

$(TARGET): $(DEP) $(OBJECTS) 
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS) $(LIBS) ${PETSC_VEC_LIB}

$(OBJECTS): %.o: %.$(EXT)
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -rf *.o $(TARGET)

%.d: %.$(EXT)
	$(CXX) $(CXXFLAGS) $< -MM -MT $(@:.d=.o) > $@

.PHONY: cleandep
cleandep:
	rm -rf *.d

.PHONY: dep
dep: $(DEP)

.PHONY: init
init:
	rm -rf bspai.o imspai.o Matrix.o
