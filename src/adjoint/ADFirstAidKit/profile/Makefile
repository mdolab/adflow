CXXFLAGS=-g -Werror -Wconversion -Wsign-promo
LDFLAGS=-lm
all: profile

build/%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $<


profile: build/reader.o build/profile.o build/string-utils.o
	$(CXX) -lm -g -o $@ $^

build/reader.o: reader.h string-utils.h
build/profile.o: reader.h

clean:
	-rm -f profile build/*
