TARGET= debugcode
# CXX = g++
CXX = mpicxx
CXXFLAGS=-Wall -O3 -fopenmp -lfftw3 -march=native
# CXXFLAGS=-Wall -Ofast -fopenmp -march=native
SRCDIR= ./src
OBJDIR= ./obj
RESTART= ./restart
SOURCE= $(wildcard $(SRCDIR)/*.cpp)
OBJECT= $(addprefix $(OBJDIR)/, $(notdir $(SOURCE:.cpp=.o)))
OUTPUT= out.$(TARGET).txt *.vtr *.vtm output stop TIEN_POD
STOP= stop

# $ export LD_LIBRARY_PATH=/home/ducnguyen/channelflow

$(TARGET): $(OBJECT)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECT) -lm

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	-mkdir -p $(OBJDIR)
	$(CXX)  $(CXXFLAGS) -o $@ -c $< 

.PHONY: clean
clean: 
	rm -rf $(TARGET) $(OUTPUT)
allclean: 
	rm -rf $(TARGET) $(OBJECT) $(OUTPUT) $(OBJDIR) $(RESTART)

run:
	make clean
	make -j4
	sh run.sh
continue:
	rm -rf $(STOP)
	make
	sh run.sh