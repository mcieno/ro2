BUILD_DIR  = $(PWD)/bin/tsp
OBJS = src/main.o src/logging.o src/tsp_fileparser.o src/tsp_solvers.o src/tsp.o src/tspplot.o          
HEADERS = ./include/*.h 
EXE = $(BUILD_DIR)
all: $(EXE) 
setting = -1   
OS := $(shell uname)


setting = 1
CPLEX_HOME = /opt/ibm/ILOG/CPLEX_Studio1210/cplex
CC = gcc
AR = ar rc
LIBS = -L${CPLEX_HOME}/lib/x86-64_linux/static_pic -L. -lcplex -lm -lpthread -ldl
INC = -I./include -I${CPLEX_HOME}/include/ilcplex


# ---------------------------------------------------------------------
# Rules
# ---------------------------------------------------------------------
CFLAGS = -Wall -O3
##CFLAGS = -Wall -g -O0 


.SUFFIXES:
.SUFFIXES: .o .c .cpp
.c.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@
.cpp.o :
	$(CC) $(CFLAGS) $(INC) -c $< -o $@

$(EXE): $(OBJS) $(LIBUTILS)
	$(CC) $(CFLAGS) -o $(EXE) $(OBJS) $(LIBS)

$(OBJS) : $(HEADERS)

$(LIBUTILS): $(OBJS_LIBUTILS)
	$(AR) $(LIBUTILS) $(OBJS_LIBUTILS)

$(LIBUTILS) : $(HEADERS_LIBUTILS)

clean:
	$(RM) $(OBJS)
	$(RM) $(OBJS_LIBUTILS)
	$(RM) $(LIBUTILS)
	$(RM) $(EXE) 
	
again:                                                               
	make clean
	make    
	
wow:
	@echo "                                      W O W W W W WWWWWWWWWWWWWWWWWWW"

who:
	@echo "you are user $(USER) with uname `uname` (OS = $(OS)) and you working with compiler setting $(setting)" 

