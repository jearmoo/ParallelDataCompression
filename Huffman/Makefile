SRCS = huffcode.cpp huffman_seq.cpp huffman_par.cpp ByteReader.cpp ByteWriter.cpp
OBJS = $(subst .cpp,.o,$(SRCS))

APP_NAME = huffcode

default: $(APP_NAME)

# compile for xeon phi
$(APP_NAME): CXX = icc -m64 -std=c++11
$(APP_NAME): CXXFLAGS = -O3 -openmp -offload-attribute-target=mic -DRUN_MIC

# compile for cpu
#cpu: CXX = g++ -m64 --std=c++11
#cpu: CXXFLAGS = -fopenmp -DDEBUG
cpu: CXX = g++ -m64 --std=c++11
cpu: CXXFLAGS = -fopenmp -O3

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(APP_NAME) $(OBJS)

cpu: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(APP_NAME) $(OBJS)

huffcode.o: huffcode.cpp huffman_seq.h huffman_par.h ByteReader.h ByteWriter.h

huffman_seq.o: huffman_seq.cpp huffman_seq.h ByteReader.h ByteWriter.h

huffman_par.o: huffman_par.cpp huffman_par.h ByteReader.h ByteWriter.h

ByteReader.o: ByteReader.cpp ByteReader.h

ByteWriter.o: ByteWriter.cpp ByteWriter.h

.PHONY: clean

clean:
	rm $(OBJS) $(APP_NAME)

submit:
	mkdir -p templates job_outputs
	cd templates && ../scripts/batch_generate.sh
	cd job_outputs && ../scripts/sub/submit.sh
