CC      = gcc
CFLAGS  = -O3 -march=native -Wall -Wextra -std=c11 -pthread
LDFLAGS = -lpthread -lm

BINS    = utpc_sieve utpc_verify_cli

.PHONY: all clean test

all: $(BINS)

utpc_sieve: utpc_sieve.c utpc.h
	$(CC) $(CFLAGS) -o $@ utpc_sieve.c $(LDFLAGS)
	@echo "Built: $@"

utpc_verify_cli: utpc_verify_cli.c utpc.h
	$(CC) $(CFLAGS) -o $@ utpc_verify_cli.c $(LDFLAGS)
	@echo "Built: $@"

test: utpc_sieve utpc_verify_cli
	@echo "\n=== Test 1: sieve to 1000 ==="
	./utpc_sieve --quiet 1000 4
	@echo "\n=== Test 2: verify known certificate ==="
	./utpc_verify_cli 96 32 2 17 5 0
	./utpc_verify_cli 146 72 2 37 9 1
	@echo "\n=== Test 3: Python verify ==="
	python3 utpc.py verify 96 32 2 17 5 0
	python3 utpc.py verify 146 72 2 37 9 1
	@echo "\n=== Test 4: Python generate + cluster (limit=500) ==="
	python3 utpc.py generate --limit 500 --clusters --quiet
	@echo "\nAll tests passed."

clean:
	rm -f $(BINS) utpc_certs.bin utpc_results.txt utpc_clusters.txt
