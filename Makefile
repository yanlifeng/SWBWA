# Toolchain
CC  = mpicc
CXX = mpicxx
AR  = ar

# SWBWA configuration
#
# EXEC_MODE:    single | cgs | cgs_cross
# CPE_ALLOCATOR: system | pool
# HOST_MALLOC_WRAPPER: 0 | 1
# HOST_MALLOC_STATS:   0 | 1 (requires HOST_MALLOC_WRAPPER=1)
# CPE_PROFILE:         0 | 1
# CPE_PROFILE_CG:      CG sampled by LWPF when CPE_PROFILE=1 (0..5)
# KSW_U8_MODE:         int32_16 | float16_16 | float16_32 (default: int32_16)
# KSW_I16_MODE:        scalar_8 | int32_8 (default: int32_8)
# MATESW_DUAL_FORWARD: 0 | 1 (150 bp same-PE forward KSW 16+16 path)
# USE_MPI:       0 | 1
# MPI-only options:
#   MPI_INPUT_MODE: static | dynamic
#   OUTPUT_MODE: split | single_unordered | discard
#   MPI_EXACT_READ_INDEX: 0 | 1 (exact n_processed for correctness checks)
EXEC_MODE           ?= single
CPE_ALLOCATOR       ?= system
HOST_MALLOC_WRAPPER ?= 1
HOST_MALLOC_STATS   ?= 0
CPE_PROFILE         ?= 0
CPE_PROFILE_CG      ?= $(if $(filter single,$(EXEC_MODE)),0,5)
KSW_U8_MODE         ?= int32_16
KSW_I16_MODE        ?= int32_8
MATESW_DUAL_FORWARD ?= 1
LWPF3_DIR            ?= /home/export/online1/mdt00/shisuan/swls-CFD/guoshi/ylf/lwpf3

USE_MPI              ?= 1
ifeq ($(USE_MPI),1)
MPI_INPUT_MODE       ?= dynamic
OUTPUT_MODE          ?= single_unordered
MPI_EXACT_READ_INDEX ?= 1
endif

ifeq ($(USE_MPI),1)
MPI_LINK_VARIANT := multi_static
else
MPI_LINK_VARIANT := wrapper_default
endif

VALID_EXEC_MODES     := single cgs cgs_cross
VALID_CPE_ALLOCATORS := system pool
VALID_BOOLEAN_VALUES := 0 1
VALID_KSW_U8_MODES    := int32_16 float16_16 float16_32
VALID_KSW_I16_MODES   := scalar_8 int32_8
VALID_MPI_INPUT_MODES := static dynamic
VALID_OUTPUT_MODES   := split single_unordered discard

ifeq ($(filter $(EXEC_MODE),$(VALID_EXEC_MODES)),)
$(error EXEC_MODE must be one of: $(VALID_EXEC_MODES))
endif
ifeq ($(filter $(CPE_ALLOCATOR),$(VALID_CPE_ALLOCATORS)),)
$(error CPE_ALLOCATOR must be one of: $(VALID_CPE_ALLOCATORS))
endif
ifeq ($(filter $(HOST_MALLOC_WRAPPER),$(VALID_BOOLEAN_VALUES)),)
$(error HOST_MALLOC_WRAPPER must be 0 or 1)
endif
ifeq ($(filter $(HOST_MALLOC_STATS),$(VALID_BOOLEAN_VALUES)),)
$(error HOST_MALLOC_STATS must be 0 or 1)
endif
ifeq ($(filter $(CPE_PROFILE),$(VALID_BOOLEAN_VALUES)),)
$(error CPE_PROFILE must be 0 or 1)
endif
ifeq ($(filter $(KSW_U8_MODE),$(VALID_KSW_U8_MODES)),)
$(error KSW_U8_MODE must be one of: $(VALID_KSW_U8_MODES))
endif
ifeq ($(filter $(KSW_I16_MODE),$(VALID_KSW_I16_MODES)),)
$(error KSW_I16_MODE must be one of: $(VALID_KSW_I16_MODES))
endif
ifeq ($(filter $(MATESW_DUAL_FORWARD),$(VALID_BOOLEAN_VALUES)),)
$(error MATESW_DUAL_FORWARD must be 0 or 1)
endif
ifeq ($(CPE_PROFILE),1)
ifeq ($(filter $(CPE_PROFILE_CG),0 1 2 3 4 5),)
$(error CPE_PROFILE_CG must be between 0 and 5)
endif
ifeq ($(EXEC_MODE),single)
ifneq ($(CPE_PROFILE_CG),0)
$(error EXEC_MODE=single requires CPE_PROFILE_CG=0)
endif
endif
endif
ifeq ($(HOST_MALLOC_STATS),1)
ifneq ($(HOST_MALLOC_WRAPPER),1)
$(error HOST_MALLOC_STATS=1 requires HOST_MALLOC_WRAPPER=1)
endif
endif
ifeq ($(filter $(USE_MPI),$(VALID_BOOLEAN_VALUES)),)
$(error USE_MPI must be 0 or 1)
endif
ifeq ($(USE_MPI),1)
ifeq ($(filter $(MPI_EXACT_READ_INDEX),$(VALID_BOOLEAN_VALUES)),)
$(error MPI_EXACT_READ_INDEX must be 0 or 1)
endif
ifeq ($(filter $(MPI_INPUT_MODE),$(VALID_MPI_INPUT_MODES)),)
$(error MPI_INPUT_MODE must be one of: $(VALID_MPI_INPUT_MODES))
endif
ifeq ($(filter $(OUTPUT_MODE),$(VALID_OUTPUT_MODES)),)
$(error OUTPUT_MODE must be one of: $(VALID_OUTPUT_MODES))
endif
endif

EXEC_MODE_VALUE_single    := SWBWA_EXEC_SINGLE_CG
EXEC_MODE_VALUE_cgs       := SWBWA_EXEC_CGS
EXEC_MODE_VALUE_cgs_cross := SWBWA_EXEC_CGS_CROSS
CPE_ALLOC_VALUE_system    := SWBWA_CPE_ALLOC_SYSTEM
CPE_ALLOC_VALUE_pool      := SWBWA_CPE_ALLOC_POOL
CPE_MALLOC_WRAPPER_system := 0
CPE_MALLOC_WRAPPER_pool   := 1
MPI_INPUT_MODE_VALUE_static  := SWBWA_MPI_INPUT_STATIC
MPI_INPUT_MODE_VALUE_dynamic := SWBWA_MPI_INPUT_DYNAMIC
OUTPUT_MODE_VALUE_split            := SWBWA_OUTPUT_SPLIT
OUTPUT_MODE_VALUE_single_unordered := SWBWA_OUTPUT_SINGLE_UNORDERED
OUTPUT_MODE_VALUE_discard          := SWBWA_OUTPUT_DISCARD
KSW_U8_MODE_VALUE_int32_16         := SWBWA_KSW_U8_INT32_16
KSW_U8_MODE_VALUE_float16_16       := SWBWA_KSW_U8_FLOAT16_16
KSW_U8_MODE_VALUE_float16_32       := SWBWA_KSW_U8_FLOAT16_32
KSW_I16_MODE_VALUE_scalar_8        := SWBWA_KSW_I16_SCALAR_8
KSW_I16_MODE_VALUE_int32_8         := SWBWA_KSW_I16_INT32_8

SWBWA_CPPFLAGS := \
	-DSWBWA_EXEC_MODE=$(EXEC_MODE_VALUE_$(EXEC_MODE)) \
	-DSWBWA_CPE_ALLOC_MODE=$(CPE_ALLOC_VALUE_$(CPE_ALLOCATOR)) \
	-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=$(HOST_MALLOC_WRAPPER) \
	-DSWBWA_ENABLE_HOST_MALLOC_STATS=$(HOST_MALLOC_STATS) \
	-DSWBWA_ENABLE_CPE_MALLOC_WRAPPER=$(CPE_MALLOC_WRAPPER_$(CPE_ALLOCATOR)) \
	-DSWBWA_ENABLE_CPE_PROFILE=$(CPE_PROFILE) \
	-DSWBWA_CPE_PROFILE_CG=$(CPE_PROFILE_CG) \
	-DSWBWA_KSW_U8_MODE=$(KSW_U8_MODE_VALUE_$(KSW_U8_MODE)) \
	-DSWBWA_KSW_I16_MODE=$(KSW_I16_MODE_VALUE_$(KSW_I16_MODE)) \
	-DSWBWA_ENABLE_MATESW_DUAL_FORWARD=$(MATESW_DUAL_FORWARD) \
	-DSWBWA_USE_MPI=$(USE_MPI)

ifeq ($(USE_MPI),1)
SWBWA_CPPFLAGS += \
	-DSWBWA_MPI_INPUT_MODE=$(MPI_INPUT_MODE_VALUE_$(MPI_INPUT_MODE)) \
	-DSWBWA_MPI_EXACT_READ_INDEX=$(MPI_EXACT_READ_INDEX) \
	-DSWBWA_OUTPUT_MODE=$(OUTPUT_MODE_VALUE_$(OUTPUT_MODE))
endif

# Compiler and linker options
OPTFLAGS  ?= -O2
WARNFLAGS ?= -Wall -Wno-unused-function
DBGFLAGS  ?= -g
INCLUDE_DIR := include
CONFIG_HEADER := $(INCLUDE_DIR)/swbwa_config.h
DEPFLAGS := -MMD -MP

CPPFLAGS += -include $(CONFIG_HEADER) $(SWBWA_CPPFLAGS)
DFLAGS   += -DHAVE_PTHREAD
INCLUDES += -I$(INCLUDE_DIR)
ifeq ($(CPE_PROFILE),1)
INCLUDES += -I$(LWPF3_DIR)
endif
CFLAGS   += $(WARNFLAGS) $(DBGFLAGS) $(OPTFLAGS) -D_GNU_SOURCE -D_FILE_OFFSET_BITS=64
CXXFLAGS += -std=c++11
LDFLAGS  +=

HOST_ARCH_FLAGS  := -mhost -fPIC -mieee -mftz -faddress_align=32
SLAVE_ARCH_FLAGS := -mslave -msimd -fPIC -mieee -mftz -faddress_align=64
HYBRID_FLAGS     := -mhybrid

LIBS := -Wl,-q -lm -lz -lpthread -lm_slave
ifeq ($(shell uname -s),Linux)
LIBS += -lrt
endif

# Targets and objects
PROG := SWBWA

HOST_DIR := src/host

LIB_OBJS := $(addprefix $(HOST_DIR)/, \
	utils.o kthread.o kstring.o ksw.o bwt.o bntseq.o bwa.o bwamem.o \
	bwamem_pair.o bwamem_extra.o malloc_wrap.o QSufSort.o bwt_gen.o \
	rope.o rle.o is.o bwtindex.o)

APP_OBJS := $(addprefix $(HOST_DIR)/, \
	bwashm.o bwase.o bwaseqio.o bwtgap.o bwtaln.o bamlite.o bwape.o \
	kopen.o pemerge.o maxk.o bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o \
	bwt_lite.o bwtsw2_chain.o fastmap.o bwtsw2_pair.o swbwa_mpi.o \
	swbwa_output.o swbwa_cpe_profile.o)

SLAVE_DIR     := src/slave
SLAVE_SOURCES := $(wildcard $(SLAVE_DIR)/*.c)
SLAVE_OBJECTS := $(SLAVE_SOURCES:.c=.o)
DEPFILES      := $(LIB_OBJS:.o=.d) $(APP_OBJS:.o=.d) \
	$(HOST_DIR)/main.d $(SLAVE_OBJECTS:.o=.d)

.PHONY: all clean depend print-config
.SUFFIXES:

all: $(PROG)

print-config:
	@echo "EXEC_MODE=$(EXEC_MODE)"
	@echo "CPE_ALLOCATOR=$(CPE_ALLOCATOR)"
	@echo "HOST_MALLOC_WRAPPER=$(HOST_MALLOC_WRAPPER)"
	@echo "HOST_MALLOC_STATS=$(HOST_MALLOC_STATS)"
	@echo "CPE_PROFILE=$(CPE_PROFILE)"
	@echo "KSW_U8_MODE=$(KSW_U8_MODE)"
ifeq ($(CPE_PROFILE),1)
	@echo "CPE_PROFILE_CG=$(CPE_PROFILE_CG)"
	@echo "LWPF3_DIR=$(LWPF3_DIR)"
endif
	@echo "USE_MPI=$(USE_MPI)"
ifeq ($(USE_MPI),1)
	@echo "MPI_INPUT_MODE=$(MPI_INPUT_MODE)"
	@echo "OUTPUT_MODE=$(OUTPUT_MODE)"
	@echo "MPI_EXACT_READ_INDEX=$(MPI_EXACT_READ_INDEX)"
endif
	@echo "MPI_LINK_VARIANT=$(MPI_LINK_VARIANT)"

# Compile rules
$(SLAVE_DIR)/%.o: $(SLAVE_DIR)/%.c $(CONFIG_HEADER)
	$(CC) $(SLAVE_ARCH_FLAGS) -c $(DEPFLAGS) $(CFLAGS) $(DFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

$(HOST_DIR)/%.o: $(HOST_DIR)/%.c $(CONFIG_HEADER)
	$(CC) $(HOST_ARCH_FLAGS) -c $(DEPFLAGS) $(CFLAGS) $(DFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

$(HOST_DIR)/%.o: $(HOST_DIR)/%.cpp $(CONFIG_HEADER)
	$(CXX) $(HOST_ARCH_FLAGS) -c $(DEPFLAGS) $(CFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

# Link and archive rules
$(PROG): libbwa.a $(APP_OBJS) $(HOST_DIR)/main.o $(SLAVE_OBJECTS)

ifeq ($(USE_MPI),1)
	@set -e; \
	link_cmd="$$( \
		$(CXX) -show $(HYBRID_FLAGS) $(CFLAGS) $(LDFLAGS) \
		$(APP_OBJS) $(HOST_DIR)/main.o $(SLAVE_OBJECTS) -o $@ -L. -lbwa $(LIBS) \
		| sed 's#/single_static#/multi_static#g' \
	)"; \
	case "$$link_cmd" in \
		*"/single_static"*) \
			echo "error: MPI link command still references single_static" >&2; \
			exit 1; \
			;; \
		*"/multi_static"*) \
			;; \
		*) \
			echo "error: MPI wrapper did not expose a multi_static library path" >&2; \
			echo "       check the output of: $(CXX) -show" >&2; \
			exit 1; \
			;; \
	esac; \
	echo "LINK $@ (MPI multi_static)"; \
	echo "$$link_cmd"; \
	eval "$$link_cmd"
else
	$(CXX) $(HYBRID_FLAGS) $(CFLAGS) $(LDFLAGS) $(APP_OBJS) $(HOST_DIR)/main.o $(SLAVE_OBJECTS) -o $@ -L. -lbwa $(LIBS)
endif

bwamem-lite: libbwa.a $(HOST_DIR)/example.o
	$(CC) $(CFLAGS) $(LDFLAGS) $(HOST_DIR)/example.o -o $@ -L. -lbwa $(LIBS)

libbwa.a: $(LIB_OBJS)
	$(AR) -csru $@ $(LIB_OBJS)

# Maintenance
clean:
	rm -f gmon.out a.out $(PROG) *.a $(HOST_DIR)/*.o $(HOST_DIR)/*.d \
		$(SLAVE_DIR)/*.o $(SLAVE_DIR)/*.d

depend:
	@echo "Dependency files are generated automatically during compilation."

-include $(DEPFILES)
