# Toolchain
CC  = mpicc
CXX = mpicxx
AR  = ar

# SWBWA configuration
#
# EXEC_MODE:    single | cgs | cgs_cross
# FORMAT_MODE:  host | cpe
# CPE_ALLOCATOR: system | pool
# HOST_MALLOC_WRAPPER: 0 | 1
# HOST_MALLOC_STATS:   0 | 1 (requires HOST_MALLOC_WRAPPER=1)
# USE_MPI:       0 | 1
# OUTPUT_MODE:   split | single_unordered
# MPI_READ_ID_SCAN: 0 | 1
EXEC_MODE           ?= single
FORMAT_MODE         ?= host
CPE_ALLOCATOR       ?= system
HOST_MALLOC_WRAPPER ?= 1
HOST_MALLOC_STATS   ?= 0
LWPF                 ?= 0

USE_MPI              ?= 1
OUTPUT_MODE          ?= split
MPI_READ_ID_SCAN     ?= 0

ifeq ($(USE_MPI),1)
MPI_LINK_VARIANT := multi_static
else
MPI_LINK_VARIANT := wrapper_default
endif

VALID_EXEC_MODES     := single cgs cgs_cross
VALID_FORMAT_MODES   := host cpe
VALID_CPE_ALLOCATORS := system pool
VALID_BOOLEAN_VALUES := 0 1
VALID_OUTPUT_MODES   := split single_unordered

ifeq ($(filter $(EXEC_MODE),$(VALID_EXEC_MODES)),)
$(error EXEC_MODE must be one of: $(VALID_EXEC_MODES))
endif
ifeq ($(filter $(FORMAT_MODE),$(VALID_FORMAT_MODES)),)
$(error FORMAT_MODE must be one of: $(VALID_FORMAT_MODES))
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
ifeq ($(HOST_MALLOC_STATS),1)
ifneq ($(HOST_MALLOC_WRAPPER),1)
$(error HOST_MALLOC_STATS=1 requires HOST_MALLOC_WRAPPER=1)
endif
endif
ifeq ($(filter $(LWPF),$(VALID_BOOLEAN_VALUES)),)
$(error LWPF must be 0 or 1)
endif
ifeq ($(filter $(USE_MPI),$(VALID_BOOLEAN_VALUES)),)
$(error USE_MPI must be 0 or 1)
endif
ifeq ($(filter $(MPI_READ_ID_SCAN),$(VALID_BOOLEAN_VALUES)),)
$(error MPI_READ_ID_SCAN must be 0 or 1)
endif
ifeq ($(filter $(OUTPUT_MODE),$(VALID_OUTPUT_MODES)),)
$(error OUTPUT_MODE must be one of: $(VALID_OUTPUT_MODES))
endif

EXEC_MODE_VALUE_single    := SWBWA_EXEC_SINGLE_CG
EXEC_MODE_VALUE_cgs       := SWBWA_EXEC_CGS
EXEC_MODE_VALUE_cgs_cross := SWBWA_EXEC_CGS_CROSS
FORMAT_MODE_VALUE_host    := SWBWA_FORMAT_HOST
FORMAT_MODE_VALUE_cpe     := SWBWA_FORMAT_CPE
CPE_ALLOC_VALUE_system    := SWBWA_CPE_ALLOC_SYSTEM
CPE_ALLOC_VALUE_pool      := SWBWA_CPE_ALLOC_POOL
CPE_MALLOC_WRAPPER_system := 0
CPE_MALLOC_WRAPPER_pool   := 1
OUTPUT_MODE_VALUE_split            := SWBWA_OUTPUT_SPLIT
OUTPUT_MODE_VALUE_single_unordered := SWBWA_OUTPUT_SINGLE_UNORDERED

SWBWA_CPPFLAGS := \
	-DSWBWA_EXEC_MODE=$(EXEC_MODE_VALUE_$(EXEC_MODE)) \
	-DSWBWA_FORMAT_MODE=$(FORMAT_MODE_VALUE_$(FORMAT_MODE)) \
	-DSWBWA_CPE_ALLOC_MODE=$(CPE_ALLOC_VALUE_$(CPE_ALLOCATOR)) \
	-DSWBWA_ENABLE_HOST_MALLOC_WRAPPER=$(HOST_MALLOC_WRAPPER) \
	-DSWBWA_ENABLE_HOST_MALLOC_STATS=$(HOST_MALLOC_STATS) \
	-DSWBWA_ENABLE_CPE_MALLOC_WRAPPER=$(CPE_MALLOC_WRAPPER_$(CPE_ALLOCATOR)) \
	-DSWBWA_ENABLE_LWPF=$(LWPF) \
	-DSWBWA_USE_MPI=$(USE_MPI) \
	-DSWBWA_ENABLE_MPI_READ_ID_SCAN=$(MPI_READ_ID_SCAN) \
	-DSWBWA_OUTPUT_MODE=$(OUTPUT_MODE_VALUE_$(OUTPUT_MODE))

# Compiler and linker options
LWPF3_DIR ?= /home/export/online1/mdt00/shisuan/sweq/ylf/someGit/lwpf3

OPTFLAGS  ?= -O2
WARNFLAGS ?= -Wall -Wno-unused-function
DBGFLAGS  ?= -g

CPPFLAGS += -include swbwa_config.h $(SWBWA_CPPFLAGS)
INCLUDES += -I$(LWPF3_DIR)
DFLAGS   += -DHAVE_PTHREAD
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

LIB_OBJS := \
	utils.o kthread.o kstring.o ksw.o bwt.o bntseq.o bwa.o bwamem.o \
	bwamem_pair.o bwamem_extra.o malloc_wrap.o QSufSort.o bwt_gen.o \
	rope.o rle.o is.o bwtindex.o

APP_OBJS := \
	bwashm.o bwase.o bwaseqio.o bwtgap.o bwtaln.o bamlite.o bwape.o \
	kopen.o pemerge.o maxk.o bwtsw2_core.o bwtsw2_main.o bwtsw2_aux.o \
	bwt_lite.o bwtsw2_chain.o fastmap.o bwtsw2_pair.o swbwa_mpi.o \
	swbwa_output.o

SLAVE_DIR     := slave
SLAVE_SOURCES := $(wildcard $(SLAVE_DIR)/*.c)
SLAVE_OBJECTS := $(SLAVE_SOURCES:.c=.o)

.PHONY: all clean depend print-config
.SUFFIXES:

all: $(PROG)

print-config:
	@echo "EXEC_MODE=$(EXEC_MODE)"
	@echo "FORMAT_MODE=$(FORMAT_MODE)"
	@echo "CPE_ALLOCATOR=$(CPE_ALLOCATOR)"
	@echo "HOST_MALLOC_WRAPPER=$(HOST_MALLOC_WRAPPER)"
	@echo "HOST_MALLOC_STATS=$(HOST_MALLOC_STATS)"
	@echo "LWPF=$(LWPF)"
	@echo "USE_MPI=$(USE_MPI)"
	@echo "OUTPUT_MODE=$(OUTPUT_MODE)"
	@echo "MPI_READ_ID_SCAN=$(MPI_READ_ID_SCAN)"
	@echo "MPI_LINK_VARIANT=$(MPI_LINK_VARIANT)"

# Compile rules
$(SLAVE_DIR)/%.o: $(SLAVE_DIR)/%.c swbwa_config.h
	$(CC) $(SLAVE_ARCH_FLAGS) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

%.o: %.c swbwa_config.h
	$(CC) $(HOST_ARCH_FLAGS) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

%.o: %.cpp swbwa_config.h
	$(CXX) $(HOST_ARCH_FLAGS) -c $(CFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

# Link and archive rules
$(PROG): libbwa.a $(APP_OBJS) main.o $(SLAVE_OBJECTS)

ifeq ($(USE_MPI),1)
	@set -e; \
	link_cmd="$$( \
		$(CXX) -show $(HYBRID_FLAGS) $(CFLAGS) $(LDFLAGS) \
		$(APP_OBJS) main.o $(SLAVE_OBJECTS) -o $@ -L. -lbwa $(LIBS) \
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
	$(CXX) $(HYBRID_FLAGS) $(CFLAGS) $(LDFLAGS) $(APP_OBJS) main.o $(SLAVE_OBJECTS) -o $@ -L. -lbwa $(LIBS)
endif

bwamem-lite: libbwa.a example.o
	$(CC) $(CFLAGS) $(LDFLAGS) example.o -o $@ -L. -lbwa $(LIBS)

libbwa.a: $(LIB_OBJS)
	$(AR) -csru $@ $(LIB_OBJS)

# Maintenance
clean:
	rm -f gmon.out *.o a.out $(PROG) *~ *.a $(SLAVE_DIR)/*.o

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) $(INCLUDES) $(CPPFLAGS) -- *.c )

# DO NOT DELETE THIS LINE -- make depend depends on it.

QSufSort.o: QSufSort.h
bamlite.o: bamlite.h malloc_wrap.h
bntseq.o: bntseq.h utils.h kseq.h malloc_wrap.h khash.h
bwa.o: bntseq.h bwa.h bwt.h ksw.h utils.h kstring.h malloc_wrap.h kvec.h
bwa.o: kseq.h
bwamem.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h ksw.h kvec.h
bwamem.o: ksort.h utils.h kbtree.h swbwa_config.h swbwa_cpe.h swbwa_cpe_layout.h
bwamem.o: swbwa_runtime.h
bwamem_extra.o: bwa.h bntseq.h bwt.h bwamem.h kstring.h malloc_wrap.h
bwamem_pair.o: kstring.h malloc_wrap.h bwamem.h bwt.h bntseq.h bwa.h kvec.h
bwamem_pair.o: utils.h ksw.h
bwape.o: bwtaln.h bwt.h kvec.h malloc_wrap.h bntseq.h utils.h bwase.h bwa.h
bwape.o: ksw.h khash.h
bwase.o: bwase.h bntseq.h bwt.h bwtaln.h utils.h kstring.h malloc_wrap.h
bwase.o: bwa.h ksw.h
bwaseqio.o: bwtaln.h bwt.h utils.h bamlite.h malloc_wrap.h kseq.h
bwashm.o: bwa.h bntseq.h bwt.h
bwt.o: utils.h bwt.h kvec.h malloc_wrap.h
bwt_gen.o: QSufSort.h malloc_wrap.h
bwt_lite.o: bwt_lite.h malloc_wrap.h
bwtaln.o: bwtaln.h bwt.h bwtgap.h utils.h bwa.h bntseq.h malloc_wrap.h
bwtgap.o: bwtgap.h bwt.h bwtaln.h malloc_wrap.h
bwtindex.o: bntseq.h bwa.h bwt.h utils.h rle.h rope.h malloc_wrap.h
bwtsw2_aux.o: bntseq.h bwt_lite.h utils.h bwtsw2.h bwt.h kstring.h
bwtsw2_aux.o: malloc_wrap.h bwa.h ksw.h kseq.h ksort.h
bwtsw2_chain.o: bwtsw2.h bntseq.h bwt_lite.h bwt.h malloc_wrap.h ksort.h
bwtsw2_core.o: bwt_lite.h bwtsw2.h bntseq.h bwt.h kvec.h malloc_wrap.h
bwtsw2_core.o: khash.h ksort.h
bwtsw2_main.o: bwt.h bwtsw2.h bntseq.h bwt_lite.h utils.h bwa.h
bwtsw2_pair.o: utils.h bwt.h bntseq.h bwtsw2.h bwt_lite.h kstring.h
bwtsw2_pair.o: malloc_wrap.h ksw.h
example.o: bwamem.h bwt.h bntseq.h bwa.h kseq.h malloc_wrap.h
fastmap.o: bwa.h bntseq.h bwt.h bwamem.h kvec.h malloc_wrap.h utils.h kseq.h
fastmap.o: swbwa_config.h swbwa_cpe.h swbwa_runtime.h
is.o: malloc_wrap.h
kopen.o: malloc_wrap.h
kstring.o: kstring.h malloc_wrap.h
ksw.o: ksw.h neon_sse.h scalar_sse.h malloc_wrap.h
main.o: kstring.h malloc_wrap.h utils.h
malloc_wrap.o: malloc_wrap.h
maxk.o: bwa.h bntseq.h bwt.h bwamem.h kseq.h malloc_wrap.h
pemerge.o: ksw.h kseq.h malloc_wrap.h kstring.h bwa.h bntseq.h bwt.h utils.h
rle.o: rle.h
rope.o: rle.h rope.h
utils.o: utils.h ksort.h malloc_wrap.h kseq.h
fastmap.o: swbwa_mpi.h swbwa_output.h
main.o: swbwa_mpi.h
swbwa_mpi.o: swbwa_config.h swbwa_mpi.h
swbwa_output.o: swbwa_config.h swbwa_mpi.h swbwa_output.h
utils.o: swbwa_mpi.h
