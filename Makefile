CC = gcc
CFLAGS = -g -march=native -mpopcnt
FAST = -O9 -march=native -mpopcnt -mlzcnt

# all: rank select
all: spider spider1L ni-spider ni-spider2L

spider: spider.c helper.h speedtest.h
	$(CC) $(FAST) spider.c helper.h speedtest.h -o spider

spider1L: spider-1L-select.c helper.h speedtest.h
	$(CC) $(FAST) spider-1L-select.c helper.h speedtest.h -o spider1L

ni-spider: ni-spider.c ni-helper.h speedtest.h
	$(CC) $(FAST) ni-spider.c ni-helper.h speedtest.h -o ni-spider

ni-spider2L: ni-spider-2L-select.c speedtest.h
	$(CC) $(FAST) ni-spider-2L-select.c speedtest.h -o ni-spider2L

clean:
	rm spider spider1L ni-spider ni-spider2L