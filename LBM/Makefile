# Include files
SOURCES=initLB.c visualLB.c boundary.c streaming.c computeCellValues.c main.c helper.c

# Compiler
# --------
CC=gcc

CFLAGS=-pedantic -Wall -std=gnu99 -fopenmp -Ofast -msse3 -mavx  -funroll-loops -g



# Linker flags
# ------------
LDFLAGS= -lpapi -fopenmp 

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=lbsim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)


$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
