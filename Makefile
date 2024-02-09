CC = gcc
CFLAGS = -g -march=native -mpopcnt
FAST = -O9 -march=native -mpopcnt -mlzcnt

# all: rank select
all: spider spider64 spider-non64 spider-non16

spider: spider.c helper.h
	$(CC) $(FAST) spider.c helper.h -o spider

spider64: spider-64-select.c helper.h
	$(CC) $(FAST) spider-64-select.c helper.h -o spider64

spider-exact: spider-exact.c
	$(CC) $(FAST) spider-exact.c -o spider-exact


debug: spider.c helper.h
	$(CC) $(CFLAGS) spider.c helper.h -o spider-debug

spider-non64: spider-non64.c
	$(CC) $(FAST) spider-non64.c -o spider-non64

spider-non16: spider-non16.c
	$(CC) $(FAST) spider-non16.c -o spider-non16

clean:
	rm spider spider64 spider-non64 spider-non16