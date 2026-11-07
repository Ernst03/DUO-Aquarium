CC ?= cc
CFLAGS ?= -O2 -std=c11 -Wall -Wextra -Wpedantic
LDLIBS = -lm

all: aquarium aquarium-gui

aquarium: aquarium.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

aquarium-gui: aquarium.c
	$(CC) $(CFLAGS) -DWITH_X11 -o $@ $< $(LDLIBS) -lX11

clean:
	rm -f aquarium aquarium-gui aquarium_events.csv aquarium_samples.csv

.PHONY: all clean
