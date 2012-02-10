CC= g++
CFLAGS=-c -Wall -ansi -pedantic-errors -Werror -g
LDFLAGS=
SOURCES=main.cpp psi.cpp theta.cpp  matrix.cpp bcg.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=main

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)
