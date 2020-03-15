BUILD_DIR  = $(PWD)/bin
SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

COSCEDIR   = /opt/ibm/ILOG/CPLEX_Studio_Community129
CPLEXDIR   = $(COSCEDIR)/cplex
CONCERTDIR = $(COSCEDIR)/concert

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CLNDIRS   = -L$(CPLEXLIBDIR)
CLNFLAGS  = -lcplex -lm -lpthread -ldl

CC         = gcc
#COPT       = -m64 -fPIC -fno-strict-aliasing
COPT       = -O3 -m64 -fno-strict-aliasing
CFLAGS     = $(COPT)  -I$(CPLEXINCDIR)

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
	$(CC) $(CFLAGS) $(CLNDIRS) -o $(BUILD_DIR)/tsp_parser src/tsp_parser.c $(CLNFLAGS)

clean:
	-@echo -e "$(FANCYLOG) Cleaning all"
	-@rm -f $(BUILD_DIR)/*
	-@rmdir $(BUILD_DIR)
