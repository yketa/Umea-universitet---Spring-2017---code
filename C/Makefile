#### BUILD PARAMETERS ####

BU=$(HOME)/code-build
OB=$(BU)/objects

#### COMPILATION PARAMETERS ####

CC=gcc
CFLAGS=-W -Wall -pedantic -std=c99 -O3
LDFLAGS=-lm

ifeq ($(TEST),yes)
	EXEC=$(BU)/test
	SRC=$(filter-out main.c, $(wildcard *.c))
else
	EXEC=$(BU)/simulation
	SRC=$(filter-out test.c, $(wildcard *.c))
endif

OBJ=$(addprefix $(OB)/, $(SRC:.c=.o))

.PHONY: all clean mrproper

#### COMPILATION #####

all: dir $(EXEC)

dir:
	@mkdir -p $(BU)
	@mkdir -p $(OB)

$(EXEC): $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS)

#### DEPENDENCIES ####

$(OB)/ellipsoids.o: ellipsoids.c ellipsoids.h init.h
	$(CC) -o $(OB)/ellipsoids.o -c ellipsoids.c $(CFLAGS)

$(OB)/maths.o: maths.c maths.h
	$(CC) -o $(OB)/maths.o -c maths.c $(CFLAGS)

$(OB)/integration.o: integration.c integration.h
	$(CC) -o $(OB)/integration.o -c integration.c $(CFLAGS)

$(OB)/measure.o: measure.c measure.h init.h
	$(CC) -o $(OB)/measure.o -c measure.c $(CFLAGS)

$(OB)/init.o: init.c init.h
	$(CC) -o $(OB)/init.o -c init.c $(CFLAGS)

$(OB)/main.o: main.c init.h ellipsoids.h param.h integration.h maths.h measure.h
	$(CC) -o $(OB)/main.o -c main.c $(CFLAGS)

$(OB)/test.o: test.c init.h ellipsoids.h param.h integration.h maths.h measure.h
	$(CC) -o $(OB)/test.o -c test.c $(CFLAGS)

#### VALGRIND ####

memcheck: dir $(OBJ)
	$(CC) -g -o $(EXEC) $(OBJ) $(LDFLAGS)
	valgrind --leak-check=yes --track-origins=yes --log-file=$(BU)/memcheck.output $(EXEC)

massif: dir $(OBJ)
	$(CC) -g -o $(EXEC) $(OBJ) $(LDFLAGS)
	valgrind --tool=massif --massif-out-file=$(BU)/massif.out $(EXEC)
	ms_print $(BU)/massif.out > $(BU)/massif.output

#### SAVE ####

source:
	@tar --exclude='.tar.gz' -zcf source_$$(date +%d-%m-%Y_%H-%M-%S).tar.gz ../$(notdir $(shell pwd)) &> /dev/null

#### CLEAN ####

clean:
	@rm -rf $(OB)

mrproper: clean
	@rm -rf $(BU)