BUILD_DIR  = $(PWD)/bin
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

COSCEDIR   = /opt/ibm/ILOG/CPLEX_Studio1210
CPLEXDIR   = $(COSCEDIR)/cplex
CONCERTDIR = $(COSCEDIR)/concert

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CONCERTINCDIR = $(CONCERTDIR)/include/ilconcert
CPLEXINCDIR   = $(CPLEXDIR)/include/ilcplex
LOCAL_INCDIR  = $(PWD)/include
TSP_SRC_FILES = src/*.c

CLNDIRS   = -L$(CPLEXLIBDIR)
CLNFLAGS  =  -lm -lpthread -ldl -lcplex

CC         = gcc
#COPT       = -O3 -m64 -fPIC
COPT       = -g3 -O0 -m64 -Wall -Werror --pedantic
CFLAGS     = $(COPT)  -I$(LOCAL_INCDIR)  -I$(CPLEXINCDIR)

FANCYLOG   = \033[96m[*]\033[0m


config:
	-@echo -e "$(FANCYLOG) BUILD_DIR     = $(BUILD_DIR)"
	-@echo -e "$(FANCYLOG) CPLEXDIR      = $(CPLEXDIR)"
	-@echo -e "$(FANCYLOG) CPLEXBINDIR   = $(CPLEXBINDIR)"
	-@echo -e "$(FANCYLOG) CPLEXLIBDIR   = $(CPLEXLIBDIR)"
	-@echo -e "$(FANCYLOG) CPLEXINCDIR   = $(CPLEXINCDIR)"
	-@echo -e "$(FANCYLOG) CONCERTDIR    = $(CONCERTDIR)"
	-@echo -e "$(FANCYLOG) CONCERTLIBDIR = $(CONCERTLIBDIR)"
	-@echo -e "$(FANCYLOG) CONCERTINCDIR = $(CONCERTINCDIR)"
	-@echo -e "$(FANCYLOG) CC            = $(CC)"
	-@echo -e "$(FANCYLOG) COPT          = $(COPT)"
	-@echo -e "$(FANCYLOG) CLNFLAGS      = $(CLNFLAGS)"

all:
	-@echo -e "$(FANCYLOG) Building all"
	mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS) $(CLNDIRS) -o $(BUILD_DIR)/tsp $(TSP_SRC_FILES) $(CLNFLAGS)

clean:
	-@echo -e "$(FANCYLOG) Cleaning all"
	rm -f $(BUILD_DIR)/*
	rmdir $(BUILD_DIR)

.PHONY: all clean config
