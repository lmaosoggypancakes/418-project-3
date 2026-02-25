APP_NAME=wireroute

OBJS=wireroute.o validate.o

CXX = g++
CXXFLAGS = -Wall -O -std=c++17 -m64 -I. -fopenmp -Wno-unknown-pragmas -pg
#-fsanitize=address
# -fsanitize=thread

all: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $<

clean:
	/bin/rm -rf *~ *.o $(APP_NAME) *.class
