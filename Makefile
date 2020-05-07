BUILD_DIR     = $(PWD)/bin
SYSTEM        = x86-64_linux
LIBFORMAT     = static_pic

COSCEDIR      = /opt/ibm/ILOG/CPLEX_Studio1210
CPLEXDIR      = $(COSCEDIR)/cplex
CONCERTDIR    = $(COSCEDIR)/concert
CONCORDEDIR   = ./concorde

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(SYSTEM)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CONCERTINCDIR = $(CONCERTDIR)/include/ilconcert
CPLEXINCDIR   = $(CPLEXDIR)/include/ilcplex
LOCAL_INCDIR  = $(PWD)/include
TSP_SRC_FILES = src/*.c

CLNDIRS       = -L$(CPLEXLIBDIR) -L$(CONCORDEDIR)
CLNFLAGS      = -Wl,--start-group -lm -lpthread -ldl -lcplex -lconcorde

CC            = gcc
COPTPERF      = -O9 -m64 -fPIC -Wall -Werror --pedantic -no-pie -fno-PIE
COPTDEBG      = -g3 -O0 -m64 -Wall -Werror --pedantic -fstack-protector-all -no-pie -fno-PIE
CFLAGSPERF    = $(COPTPERF)  -I$(LOCAL_INCDIR)  -I$(CPLEXINCDIR)
CFLAGSDEBG    = $(COPTDEBG)  -I$(LOCAL_INCDIR)  -I$(CPLEXINCDIR)

FANCYLOG      = \033[96m[*]\033[0m


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
	$(CC) $(CFLAGSPERF) $(CLNDIRS) -o $(BUILD_DIR)/tsp $(TSP_SRC_FILES) $(CLNFLAGS)

debug:
	-@echo -e "$(FANCYLOG) Building all for debugging"
	mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGSDEBG) $(CLNDIRS) -o $(BUILD_DIR)/tsp $(TSP_SRC_FILES) $(CLNFLAGS)

clean:
	-@echo -e "$(FANCYLOG) Cleaning all"
	rm -f $(BUILD_DIR)/*
	rmdir $(BUILD_DIR) || exit 0

.PHONY: debug perf clean config
