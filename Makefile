# Taken from http://www.partow.net/programming/makefile/index.html
MAKEFLAGS += --silent
CXX      := -g++
CXXFLAGS := -std=c++11 -fopenmp -Isrc/ -MD

LDFLAGS  := 
BUILD    := bin
LIBFLAGS :=
OBJ_DIR  := $(BUILD)
APP_DIR  := $(BUILD)
TARGET   := cosims
INCLUDE  :=
SRC      :=                      \
	$(wildcard src/*.cpp)         \

OBJECTS := $(SRC:%.cpp=$(OBJ_DIR)/%.o)

# compile wiht little optimizaiton and no MKL libraries
ifeq ($(config), debug)
	CXX := -g++
	CXXFLAGS := $(CXXFLAGS) -O1
	TARGET := $(TARGET)_debug

# use Intel MKL libraries and Intel Compiler
else
	CXX := -g++
	CXXFLAGS := $(CXXFLAGS) -O3 -static
	TARGET := $(TARGET)
endif



all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o: %.cpp
	@echo "Building " $<
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $< $(LIBFLAGS) $(INTELLINK)

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@echo "\nLinking " $(TARGET)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o $(APP_DIR)/$(TARGET) $(OBJECTS) $(LIBFLAGS)
	@echo "Done"

.PHONY: all build clean debug release

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

debug: CXXFLAGS += -DDEBUG -g
debug: all

release: CXXFLAGS += -O2
release: all

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

-include $(OBJECTS:.o=.d)
