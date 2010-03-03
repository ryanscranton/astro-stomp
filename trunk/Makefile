CC = g++
LD = g++
CCFLAGS = -O6 -Wall
LDFLAGS = -lm -L./ -lstomp -lgflags
LIB=libstomp.a

CORE = stomp_core
ANGULAR_COORDINATE = stomp_angular_coordinate
ANGULAR_BIN = stomp_angular_bin
ANGULAR_CORRELATION = stomp_angular_correlation
PIXEL = stomp_pixel
SCALAR_PIXEL = stomp_scalar_pixel
TREE_PIXEL = stomp_tree_pixel
BASE_MAP = stomp_base_map
MAP = stomp_map
SCALAR_MAP = stomp_scalar_map
TREE_MAP = stomp_tree_map
FOOTPRINT = stomp_footprint
UTIL = stomp_util

TEST_EXEC = stomp_unit_test
GENRAND_EXEC = stomp_genrand

default: lib

all: lib test genrand

lib: $(LIB)

$(LIB): $(CORE).o $(ANGULAR_COORDINATE).o $(ANGULAR_BIN).o $(ANGULAR_CORRELATION).o $(PIXEL).o $(SCALAR_PIXEL).o $(TREE_PIXEL).o $(BASE_MAP).o $(MAP).o $(SCALAR_MAP).o $(TREE_MAP).o $(FOOTPRINT).o $(UTIL).o
	@ echo linking $(LIB)
	ar -rv $(LIB) $(CORE).o $(ANGULAR_COORDINATE).o $(ANGULAR_BIN).o $(ANGULAR_CORRELATION).o $(PIXEL).o $(SCALAR_PIXEL).o $(TREE_PIXEL).o $(BASE_MAP).o $(MAP).o $(SCALAR_MAP).o $(TREE_MAP).o $(FOOTPRINT).o $(UTIL).o

$(CORE).o: $(CORE).cc $(CORE).h
	$(CC) -c $(CCFLAGS) $(CORE).cc

$(ANGULAR_COORDINATE).o: $(ANGULAR_COORDINATE).cc $(ANGULAR_COORDINATE).h $(CORE).h $(PIXEL).h
	$(CC) -c $(CCFLAGS) $(ANGULAR_COORDINATE).cc

$(ANGULAR_BIN).o: $(ANGULAR_BIN).cc $(ANGULAR_BIN).h $(CORE).h
	$(CC) -c $(CCFLAGS) $(ANGULAR_BIN).cc

$(ANGULAR_CORRELATION).o: $(ANGULAR_CORRELATION).cc $(ANGULAR_CORRELATION).h $(CORE).h $(ANGULAR_COORDINATE).h $(ANGULAR_BIN).h $(MAP).h $(SCALAR_MAP).h $(TREE_MAP).h
	$(CC) -c $(CCFLAGS) $(ANGULAR_CORRELATION).cc

$(PIXEL).o: $(PIXEL).cc $(PIXEL).h $(CORE).h $(ANGULAR_COORDINATE).h MersenneTwister.h $(ANGULAR_BIN).h
	$(CC) -c $(CCFLAGS) $(PIXEL).cc

$(SCALAR_PIXEL).o: $(SCALAR_PIXEL).cc $(SCALAR_PIXEL).h $(CORE).h $(PIXEL).h
	$(CC) -c $(CCFLAGS) $(SCALAR_PIXEL).cc

$(TREE_PIXEL).o: $(TREE_PIXEL).cc $(TREE_PIXEL).h $(CORE).h $(PIXEL).h $(ANGULAR_COORDINATE).h $(ANGULAR_BIN).h $(ANGULAR_CORRELATION).h
	$(CC) -c $(CCFLAGS) $(TREE_PIXEL).cc

$(BASE_MAP).o: $(BASE_MAP).cc $(BASE_MAP).h $(CORE).h $(PIXEL).h
	$(CC) -c $(CCFLAGS) $(BASE_MAP).cc

$(MAP).o: $(MAP).cc $(MAP).h $(CORE).h $(ANGULAR_COORDINATE).h $(PIXEL).h $(BASE_MAP).h
	$(CC) -c $(CCFLAGS) $(MAP).cc

$(SCALAR_MAP).o: $(SCALAR_MAP).cc $(SCALAR_MAP).h $(CORE).h $(SCALAR_PIXEL).h $(ANGULAR_BIN).h $(BASE_MAP).h
	$(CC) -c $(CCFLAGS) $(SCALAR_MAP).cc

$(TREE_MAP).o: $(TREE_MAP).cc $(TREE_MAP).h $(CORE).h $(TREE_PIXEL).h $(BASE_MAP).h $(ANGULAR_COORDINATE).h $(ANGULAR_BIN).h $(ANGULAR_CORRELATION).h
	$(CC) -c $(CCFLAGS) $(TREE_MAP).cc

$(FOOTPRINT).o: $(FOOTPRINT).cc $(FOOTPRINT).h $(CORE).h $(PIXEL).h $(MAP).h $(ANGULAR_COORDINATE).h
	$(CC) -c $(CCFLAGS) $(FOOTPRINT).cc

$(UTIL).o: $(UTIL).cc $(UTIL).h $(CORE).h
	$(CC) -c $(CCFLAGS) $(UTIL).cc

# unit tests
test: $(TEST_EXEC).o $(CORE)_test.o $(ANGULAR_COORDINATE)_test.o $(ANGULAR_CORRELATION)_test.o $(PIXEL)_test.o $(SCALAR_PIXEL)_test.o $(TREE_PIXEL)_test.o $(MAP)_test.o $(SCALAR_MAP)_test.o $(TREE_MAP)_test.o $(FOOTPRINT)_test.o $(UTIL)_test.o
	$(LD) $(TEST_EXEC).o $(CORE)_test.o $(ANGULAR_COORDINATE)_test.o $(ANGULAR_CORRELATION)_test.o $(PIXEL)_test.o $(SCALAR_PIXEL)_test.o $(TREE_PIXEL)_test.o $(MAP)_test.o $(SCALAR_MAP)_test.o $(TREE_MAP)_test.o $(FOOTPRINT)_test.o $(UTIL)_test.o $(LDFLAGS) -o $(TEST_EXEC)

$(TEST_EXEC).o: $(LIB) $(TEST_EXEC).cc
	$(CC) -c $(CCFLAGS) $(TEST_EXEC).cc

$(CORE)_test.o: $(LIB) $(CORE)_test.cc $(CORE)_test.h
	$(CC) -c $(CCFLAGS) $(CORE)_test.cc

$(ANGULAR_COORDINATE)_test.o: $(ANGULAR_COORDINATE)_test.cc $(ANGULAR_COORDINATE)_test.h $(LIB)
	$(CC) -c $(CCFLAGS) $(ANGULAR_COORDINATE)_test.cc

$(ANGULAR_CORRELATION)_test.o: $(LIB) $(ANGULAR_CORRELATION)_test.cc $(ANGULAR_CORRELATION)_test.h
	$(CC) -c $(CCFLAGS) $(ANGULAR_CORRELATION)_test.cc

$(PIXEL)_test.o: $(LIB) $(PIXEL)_test.cc $(PIXEL)_test.h
	$(CC) -c $(CCFLAGS) $(PIXEL)_test.cc

$(SCALAR_PIXEL)_test.o: $(LIB) $(SCALAR_PIXEL)_test.cc $(SCALAR_PIXEL)_test.h
	$(CC) -c $(CCFLAGS) $(SCALAR_PIXEL)_test.cc

$(TREE_PIXEL)_test.o: $(LIB) $(TREE_PIXEL)_test.cc $(TREE_PIXEL)_test.h
	$(CC) -c $(CCFLAGS) $(TREE_PIXEL)_test.cc

$(MAP)_test.o: $(LIB) $(MAP)_test.cc $(MAP)_test.cc
	$(CC) -c $(CCFLAGS) $(MAP)_test.cc

$(SCALAR_MAP)_test.o: $(LIB) $(SCALAR_MAP)_test.cc $(SCALAR_MAP)_test.h
	$(CC) -c $(CCFLAGS) $(SCALAR_MAP)_test.cc

$(TREE_MAP)_test.o: $(LIB) $(TREE_MAP)_test.cc $(TREE_MAP)_test.h
	$(CC) -c $(CCFLAGS) $(TREE_MAP)_test.cc

$(FOOTPRINT)_test.o: $(LIB) $(FOOTPRINT)_test.cc $(FOOTPRINT)_test.h 
	$(CC) -c $(CCFLAGS) $(FOOTPRINT).cc

$(UTIL)_test.o: $(LIB) $(UTIL)_test.cc $(UTIL)_test.h
	$(CC) -c $(CCFLAGS) $(UTIL)_test.cc


# executable to generate random points and write various
# formats to stdout
genrand: $(GENRAND_EXEC).o
	$(LD) $(GENRAND_EXEC).o $(LDFLAGS) -o $(GENRAND_EXEC)

$(GENRAND_EXEC).o: $(LIB) $(GENRAND_EXEC).cc
	$(CC) -c $(CCFLAGS) $(GENRAND_EXEC).cc


clean:
	rm -f *.o *.a $(TEST_EXEC) $(GENRAND_EXEC)



