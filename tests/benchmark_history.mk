# Test-only replacement main, read by GNU make through MAKEFILES.
# Keep the regular link and build.sh's two-pass cross relocation workflow.
.DEFAULT_GOAL := all
HISTORY_BENCH_OBJECTS := tests/benchmark_history_runtime.o tests/benchmark_pread_cache.o

src/host/main.o: tests/benchmark_mpi_gzip_history.cpp $(HISTORY_BENCH_OBJECTS)
	$(CXX) $(HOST_ARCH_FLAGS) -c $(DEPFLAGS) $(CFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(CPPFLAGS) -Isrc/gzip -DSWBWA_HISTORY_APPLICATION_BACKEND=1 $< -o $@

tests/benchmark_history_runtime.o: tests/benchmark_history_runtime.c
	$(CC) $(HOST_ARCH_FLAGS) -c $(DEPFLAGS) $(CFLAGS) $(DFLAGS) $(INCLUDES) $(CPPFLAGS) $< -o $@

tests/benchmark_pread_cache.o: tests/benchmark_pread_cache.cpp
	$(CXX) $(HOST_ARCH_FLAGS) -c $(DEPFLAGS) $(CFLAGS) $(CXXFLAGS) -D_LARGEFILE64_SOURCE $< -o $@

.PHONY: history-benchmark-clean
clean: history-benchmark-clean
history-benchmark-clean:
	$(RM) $(HISTORY_BENCH_OBJECTS) $(HISTORY_BENCH_OBJECTS:.o=.d)

-include $(HISTORY_BENCH_OBJECTS:.o=.d)
