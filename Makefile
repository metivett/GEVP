CC = g++

CFLAGS = -O2 -std=c++11
LDFLAGS = -fopenmp

INCLUDES = -I./include
INCLUDES += -I/usr/local/include/
INCLUDES += -I/home/thibaut/workspace/LQCDAnalysis/include
INCLUDES += -I/home/thibaut/workspace/LQCDAnalysis/
INCLUDES += -I/home/thibaut/workspace/LQCDAnalysis/utils/include

LIBS = -L"/home/thibaut/workspace/LQCDAnalysis/" -lLQCDAnalysis
LIBS += -L"/home/thibaut/workspace/LQCDAnalysis/utils" -lutils
LIBS += -lMinuit2
LIBS += -lboost_program_options
LIBS += -lgsl -lgslcblas

SRC_DIR = ./src
OBJ_FILES = Z001.o gevp.o plateau.o utils.o analyze.o main.o

EXE_NAME = gevp

all : $(EXE_NAME)

clean :
	rm $(EXE_NAME) $(OBJ_FILES)

$(EXE_NAME): $(OBJ_FILES)
	@echo 'Building target $@'
	@echo 'Invoking GCC C++ Linker'
	$(CC) $(LDFLAGS) -o $(EXE_NAME) $(OBJ_FILES) $(LIBS)
	@echo 'Finished building target $@'
	@echo ' '

%.o: $(SRC_DIR)/%.cpp
	@echo 'Building file $<'
	@echo 'Invoking GCC C++ Compiler'
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
	@echo 'Finished building: $<'
	@echo ' '

%.o: $(SRC_DIR)/%.cc
	@echo 'Building file $<'
	@echo 'Invoking GCC C++ Compiler'
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
	@echo 'Finished building: $<'
	@echo ' '

%.o: $(SRC_DIR)/%.c
	@echo 'Building file $<'
	@echo 'Invoking GCC C++ Compiler'
	gcc $(CFLAGS) $(INCLUDES) -o $@ -c $<
	@echo 'Finished building: $<'
	@echo ' '

.PHONY: clean $(EXE_NAME)