# Builds the C helpers of tools/evaluate.py into tools/bin/ (git-ignored).
#   make -C tools/c
CC ?= cc
CFLAGS ?= -O3 -Wall -Wextra -Wno-unused-parameter
BIN := ../bin

all: $(BIN)/libvntreval.so

$(BIN)/libvntreval.so: vntreval.c
	mkdir -p $(BIN)
	$(CC) $(CFLAGS) -shared -fPIC -o $@ $<

clean:
	rm -f $(BIN)/libvntreval.so

.PHONY: all clean
