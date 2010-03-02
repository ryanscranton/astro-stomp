CC = g++
LD = g++
CCFLAGS = -O6 -Wall
LDFLAGS = -lm -L./ -lstomp -lgflags
LIB=libstomp.so

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
COSMOLOGY = stomp_cosmology
UTIL = stomp_util

TEST_EXEC = stomp_unit_test
GENRAND_EXEC = stomp_genrand

default: lib

all: lib test genrand

lib: $(LIB)

$(LIB): $(CORE).o $(ANGULAR_COORDINATE).o $(ANGULAR_BIN).o $(ANGULAR_CORRELATION).o $(PIXEL).o $(SCALAR_PIXEL).o $(TREE_PIXEL).o $(BASE_MAP).o $(MAP).o $(SCALAR_MAP).o $(TREE_MAP).o $(FOOTPRINT).o $(UTIL).o
	@ echo linking $(LIB)
	$(CC) -shared -o $(LIB) $(CORE).o $(ANGULAR_COORDINATE).o $(ANGULAR_BIN).o $(ANGULAR_CORRELATION).o $(PIXEL).o $(SCALAR_PIXEL).o $(TREE_PIXEL).o $(BASE_MAP).o $(MAP).o $(SCALAR_MAP).o $(TREE_MAP).o $(FOOTPRINT).o $(UTIL).o

$(CORE).o: $(CORE).cc $(CORE).h
	$(CC) -c -fPIC $(CCFLAGS) $(CORE).cc

$(ANGULAR_COORDINATE).o: $(ANGULAR_COORDINATE).cc $(ANGULAR_COORDINATE).h $(CORE).h $(PIXEL).h
	$(CC) -c -fPIC $(CCFLAGS) $(ANGULAR_COORDINATE).cc

$(ANGULAR_BIN).o: $(ANGULAR_BIN).cc $(ANGULAR_BIN).h $(CORE).h
	$(CC) -c -fPIC $(CCFLAGS) $(ANGULAR_BIN).cc

$(ANGULAR_CORRELATION).o: $(ANGULAR_CORRELATION).cc $(ANGULAR_CORRELATION).h $(CORE).h $(ANGULAR_COORDINATE).h $(ANGULAR_BIN).h $(MAP).h $(SCALAR_MAP).h $(TREE_MAP).h
	$(CC) -c -fPIC $(CCFLAGS) $(ANGULAR_CORRELATION).cc

$(PIXEL).o: $(PIXEL).cc $(PIXEL).h $(CORE).h $(ANGULAR_COORDINATE).h MersenneTwister.h $(ANGULAR_BIN).h
	$(CC) -c -fPIC $(CCFLAGS) $(PIXEL).cc

$(SCALAR_PIXEL).o: $(SCALAR_PIXEL).cc $(SCALAR_PIXEL).h $(CORE).h $(PIXEL).h
	$(CC) -c -fPIC $(CCFLAGS) $(SCALAR_PIXEL).cc

$(TREE_PIXEL).o: $(TREE_PIXEL).cc $(TREE_PIXEL).h $(CORE).h $(PIXEL).h $(ANGULAR_COORDINATE).h $(ANGULAR_BIN).h $(ANGULAR_CORRELATION).h
	$(CC) -c -fPIC $(CCFLAGS) $(TREE_PIXEL).cc

$(BASE_MAP).o: $(BASE_MAP).cc $(BASE_MAP).h $(CORE).h $(PIXEL).h
	$(CC) -c -fPIC $(CCFLAGS) $(BASE_MAP).cc

$(MAP).o: $(MAP).cc $(MAP).h $(CORE).h $(ANGULAR_COORDINATE).h $(PIXEL).h $(BASE_MAP).h
	$(CC) -c -fPIC $(CCFLAGS) $(MAP).cc

$(SCALAR_MAP).o: $(SCALAR_MAP).cc $(SCALAR_MAP).h $(CORE).h $(SCALAR_PIXEL).h $(ANGULAR_BIN).h $(BASE_MAP).h
	$(CC) -c -fPIC $(CCFLAGS) $(SCALAR_MAP).cc

$(TREE_MAP).o: $(TREE_MAP).cc $(TREE_MAP).h $(CORE).h $(TREE_PIXEL).h $(BASE_MAP).h $(ANGULAR_COORDINATE).h $(ANGULAR_BIN).h $(ANGULAR_CORRELATION).h
	$(CC) -c -fPIC $(CCFLAGS) $(TREE_MAP).cc

$(FOOTPRINT).o: $(FOOTPRINT).cc $(FOOTPRINT).h $(CORE).h $(PIXEL).h $(MAP).h $(ANGULAR_COORDINATE).h
	$(CC) -c -fPIC $(CCFLAGS) $(FOOTPRINT).cc

$(UTIL).o: $(UTIL).cc $(UTIL).h $(CORE).h
	$(CC) -c -fPIC $(CCFLAGS) $(UTIL).cc

# unit tests
test: $(TEST_EXEC).o
	$(LD) $(TEST_EXEC).o $(LDFLAGS) -o $(TEST_EXEC)

$(TEST_EXEC).o: $(LIB) $(TEST_EXEC).cc
	$(CC) -c -fPIC $(CCFLAGS) $(TEST_EXEC).cc


# executable to generate random points and write various
# formats to stdout
genrand: $(GENRAND_EXEC).o
	$(LD) $(GENRAND_EXEC).o $(LDFLAGS) -o $(GENRAND_EXEC)

$(GENRAND_EXEC).o: $(LIB) $(GENRAND_EXEC).cc
	$(CC) -c $(CCFLAGS) $(GENRAND_EXEC).cc


clean:
	rm -f *.o *.a $(TEST_EXEC) $(GENRAND_EXEC)



