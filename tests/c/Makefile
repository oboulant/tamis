CC=gcc
CFLAGS=-Wall
INC_EBCD_FLAGS := -I../../src/tamis/_tamis/
INC_UNITY_FLAGS := -Iunity/src
OBJ = ./ebcd/ebcd.o
LIBRARIES = -lm

%.o: %.c %.h
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJ_TEST)
	$(CC) -c -o ./TestEbcd.o ./TestEbcd.c $(INC_UNITY_FLAGS) $(INC_EBCD_FLAGS) $(CFLAGS)
	$(CC) -c -o ./unity/src/unity.o ./unity/src/unity.c $(INC_UNITY_FLAGS) $(CFLAGS)
	$(CC) -o ./test ./TestEbcd.o ./unity/src/unity.o $(CFLAGS) $(LIBRARIES)

.PHONY: clean

clean:
	find . -type f -name '*.o' -exec rm {} +
	rm -f test
