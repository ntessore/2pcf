-include local.mk

CFLAGS += -std=c99 -Wall -Wextra -Werror -Wno-unknown-pragmas -pedantic
LDFLAGS += 
LDLIBS += -lm

## system-specific libraries
#LDLIBS_Windows_NT += (BLAS lib?)
LDLIBS_Linux += -llapack -lblas
LDLIBS_Darwin += -framework Accelerate

ifdef DEBUG
CFLAGS += -O0 -g -DDEBUG
else
CFLAGS += -Ofast -march=native -mfpmath=sse
endif

ifdef OPENMP
CFLAGS += -fopenmp
endif

ifndef OS
OS = $(shell uname -s)
endif

LDLIBS += $(LDLIBS_$(OS))

.PHONY: all clean

all: 2pcf

clean:
	$(RM) 2pcf

2pcf: 2pcf.c io.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $< $(LDLIBS)
