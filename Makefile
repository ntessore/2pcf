-include local.mk

CFLAGS += -std=c99 -Wall -Wextra -Werror -Wno-unknown-pragmas -pedantic
LDFLAGS += 
LDLIBS += -lm

ifdef DEBUG
CFLAGS += -O0 -g -DDEBUG
else
CFLAGS += -Ofast
endif

.PHONY: all clean

all: 2pcf

clean:
	$(RM) 2pcf

2pcf: 2pcf.c
	$(CC) $(CFLAGS) $(LDFLAGS) $(CPPFLAGS) -o $@ $^ $(LDLIBS)
