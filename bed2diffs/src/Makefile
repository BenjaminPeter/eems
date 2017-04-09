
.PHONY: all

PLINKIO = /home/bpeter/programs/libplinkio/src/
EXE = bed2diffs_v1
OBJ = bed2diffs_v1.o data.o

CXXFLAGS = -I${PLINKIO}/include -O3 -Wall -Werror -fopenmp -static
LDFLAGS = -lplinkio -L${PLINKIO}/lib 

all: $(EXE)
$(EXE): $(OBJ)
	$(CXX)  $(CXXFLAGS) -o $(EXE) $^ $(LDFLAGS)
$(OBJ): %.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
clean:
	rm -f $(OBJ) $(EXE) *~
