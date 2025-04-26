# Compiler and flags
CC := gcc
CXX := g++
# I do not add std flags because if we add libraries that are not compatible, it will break the compilation
CXXFLAGS := -Wall -Wextra -Wpedantic -Wshadow -Wformat=2 -Wcast-align -Wconversion -Wsign-conversion -Wnull-dereference -g3 -Ofast 
CFLAGS := -Wall -Wextra -pedantic -Ofast

# Libraries
LIBS := -lm -lmatplot -lpthread

# Folders
SRC := src
INCLUDE := -Iinclude -Ilibs
BIN := bin
LIB_DIR := libs

# Executable name
TARGET := main

# Sources and objects
SOURCES := $(wildcard $(SRC)/*.c) $(wildcard $(SRC)/*.cpp)
SOURCES_LIB := $(wildcard $(LIB_DIR)/*.c) $(wildcard $(LIB_DIR)/*.cpp)
OBJECTS := $(patsubst $(SRC)/%.cpp, $(BIN)/%.o, $(filter %.cpp, $(SOURCES))) \
           $(patsubst $(SRC)/%.c, $(BIN)/%.o, $(filter %.c, $(SOURCES)))

# Library objects (will be stored separately)
LIB_OBJECTS := $(patsubst $(LIB_DIR)/%.c, $(LIB_DIR)/%.o, $(filter %.c, $(SOURCES_LIB))) \
               $(patsubst $(LIB_DIR)/%.cpp, $(LIB_DIR)/%.o, $(filter %.cpp, $(SOURCES_LIB)))

# Headers
HEADERS := $(wildcard $(INCLUDE)/*.h) $(wildcard $(INCLUDE)/*.hpp) \
           $(wildcard $(LIB_DIR)/*.h) $(wildcard $(LIB_DIR)/*.hpp)

# .PHONY target specifies that all and clean are not real files, but are just targets that don't produce output files.
.PHONY: all clean build-libs

# Default target
all: $(BIN)/$(TARGET)

# Build library objects first
build-libs: $(LIB_OBJECTS)

$(BIN)/$(TARGET): $(OBJECTS) $(LIB_OBJECTS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

$(BIN)/%.o: $(SRC)/%.cpp $(HEADERS) | $(BIN)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

$(BIN)/%.o: $(SRC)/%.c $(HEADERS) | $(BIN)
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Library objects compilation rules
$(LIB_DIR)/%.o: $(LIB_DIR)/%.cpp $(HEADERS) | $(LIB_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

$(LIB_DIR)/%.o: $(LIB_DIR)/%.c $(HEADERS) | $(LIB_DIR)
	$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

# Ensure directories exist
$(BIN):
	@mkdir -p $(BIN)

# Clean target (doesn't remove library objects)
clean:
	@rm -rf $(BIN)/*.o $(BIN)/$(TARGET)

run: all
	@./$(BIN)/$(TARGET); python -m src.main
