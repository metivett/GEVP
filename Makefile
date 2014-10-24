# CC = g++-4.9
CC = g++-4.7
# CC = clang

CFLAGS = -O2 -std=c++11 -g
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
LIBS += -lboost_system
LIBS += -lboost_filesystem
LIBS += -lgsl -lgslcblas
LIBS += -lgrace_np
LIBS += -lm

SRC_DIR = ./src
COMMON_OBJ_FILES = plateau.o utils.o
GEVP_OBJ_FILES = Z001.o gevp.o analyze.o gevp_main.o
SET_SCALE_OBJ_FILES = set_scale.o set_scale_main.o
EXTRACT_OBJ_FILES = extract.o extract_main.o

.PHONY: all clean

all: gevp set_scale extract

clean:
	rm -f gevp set_scale extract
	rm -f $(COMMON_OBJ_FILES)
	rm -f $(GEVP_OBJ_FILES)
	rm -f $(SET_SCALE_OBJ_FILES)
	rm -f $(EXTRACT_OBJ_FILES)

gevp: $(COMMON_OBJ_FILES) $(GEVP_OBJ_FILES)
	@echo 'Building target $@'
	@echo 'Invoking GCC C++ Linker'
	$(CC) $(LDFLAGS) -o $@ $(COMMON_OBJ_FILES) $(GEVP_OBJ_FILES) $(LIBS)
	@echo 'Finished building target $@'
	@echo ' '

set_scale: $(COMMON_OBJ_FILES) $(SET_SCALE_OBJ_FILES)
	@echo 'Building target $@'
	@echo 'Invoking GCC C++ Linker'
	$(CC) $(LDFLAGS) -o $@ $(COMMON_OBJ_FILES) $(SET_SCALE_OBJ_FILES) $(LIBS)
	@echo 'Finished building target $@'
	@echo ' '

extract: $(COMMON_OBJ_FILES) $(EXTRACT_OBJ_FILES)
	@echo 'Building target $@'
	@echo 'Invoking GCC C++ Linker'
	$(CC) $(LDFLAGS) -o $@ $(COMMON_OBJ_FILES) $(EXTRACT_OBJ_FILES) $(LIBS)
	@echo 'Finished building target $@'
	@echo ' '

$(COMMON_OBJ_FILES): %.o: $(SRC_DIR)/common/%.cpp
	@echo 'Building file $<'
	@echo 'Invoking GCC C++ Compiler'
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
	@echo 'Finished building: $<'
	@echo ' '

$(GEVP_OBJ_FILES): %.o: $(SRC_DIR)/gevp/%.cpp
	@echo 'Building file $<'
	@echo 'Invoking GCC C++ Compiler'
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
	@echo 'Finished building: $<'
	@echo ' '

$(SET_SCALE_OBJ_FILES): %.o: $(SRC_DIR)/set_scale/%.cpp
	@echo 'Building file $<'
	@echo 'Invoking GCC C++ Compiler'
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
	@echo 'Finished building: $<'
	@echo ' '

$(EXTRACT_OBJ_FILES): %.o: $(SRC_DIR)/extract/%.cpp
	@echo 'Building file $<'
	@echo 'Invoking GCC C++ Compiler'
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ -c $<
	@echo 'Finished building: $<'
	@echo ' '