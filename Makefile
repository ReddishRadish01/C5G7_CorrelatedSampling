NVCC ?= nvcc
CXXSTD ?= c++17
CUDA_ARCH ?= native
BUILD_DIR ?= build
TARGET ?= $(BUILD_DIR)/c5g7
THREADS_PER_BLOCK ?= 32
DEVICE_LTO ?= 1
FAST_MATH ?= 0

SRCS := main.cu Core.cu Cycle.cu Neutron.cu RNG.cu Tally.cu XSParser.cu
HDRS := $(wildcard *.cuh)

NVCCFLAGS ?= -O3 -DNDEBUG -std=$(CXXSTD) -rdc=true -arch=$(CUDA_ARCH) \
	-DC5G7_THREADS_PER_BLOCK=$(THREADS_PER_BLOCK) -split-compile=0 \
	-Xptxas -O3 -Xcompiler=-O3

ifeq ($(DEVICE_LTO),1)
NVCCFLAGS += -dlto
endif

ifeq ($(FAST_MATH),1)
NVCCFLAGS += --use_fast_math
endif

.PHONY: all clean run-critical-sweep run-critical-binary-search run-critical-newton-cs

all: $(TARGET)

$(TARGET): $(SRCS) $(HDRS) | $(BUILD_DIR)
	mkdir -p $(dir $@)
	$(NVCC) $(NVCCFLAGS) $(SRCS) -o $@

$(BUILD_DIR):
	mkdir -p $@

run-critical-sweep: $(TARGET)
	bash scripts/run_critical_rod_sweep.sh

run-critical-binary-search: $(TARGET)
	bash scripts/run_critical_rod_binary_search.sh

run-critical-newton-cs: $(TARGET)
	bash scripts/run_critical_rod_newton_cs.sh

clean:
	rm -rf $(BUILD_DIR)
