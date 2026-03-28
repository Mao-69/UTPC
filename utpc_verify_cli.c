/*
 * utpc_verify_cli.c — Command-line UTPC certificate verifier
 *
 * Usage:
 *   ./utpc_verify_cli <n> <phi> <mu> <J> <dr> <sg>
 *   ./utpc_verify_cli --hex <27-byte-hexstring>
 *   ./utpc_verify_cli --file <certs.bin>     (verify all certs in file)
 *
 * Exit code: 0 if all certificates valid, 1 if any invalid.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utpc.h"

#define GREEN  "\033[92m"
#define RED    "\033[91m"
#define YELLOW "\033[93m"
#define RESET  "\033[0m"

static void print_result(const utpc_t *c, const utpc_verify_result_t *r,
                         int verbose) {
    const char *status = r->valid ? GREEN "VALID" RESET : RED "INVALID" RESET;

    printf("UTPC { n=%lu phi=%lu mu=%u J=%lu dr=%u sg=%u }  %s\n",
           c->n, c->phi, c->mu, c->J, c->dr, c->sg, status);

    if (verbose) {
        #define MARK(v) ((v) == 1 ? GREEN"✓"RESET : (v) == 0 ? RED"✗"RESET : YELLOW"−"RESET)
        printf("  [1+2] phi_correct   %s\n", MARK(r->phi_correct));
        printf("  [ 3 ] dr_correct    %s\n", MARK(r->dr_correct));
        printf("  [ 4 ] mu_correct    %s\n", MARK(r->mu_correct));
        printf("  [ 5 ] divisibility  %s\n", MARK(r->divisible));
        printf("  [ 6 ] J_formula     %s\n", MARK(r->J_correct));
        printf("  [ 7 ] J_prime       %s\n", MARK(r->J_prime));
        printf("  [ 8 ] untouchable   %s  (not checked — requires sieve)\n",
               MARK(r->n_untouchable));
        printf("  [ 9 ] sg_correct    %s\n", MARK(r->sg_correct));
        printf("\n");
        #undef MARK
    }
}

static int parse_hex(const char *hex, utpc_t *cert) {
    size_t len = strlen(hex);
    /* Allow spaces in hex string */
    char clean[UTPC_WIRE_SIZE * 2 + 1];
    int ci = 0;
    for (size_t i = 0; i < len && ci < UTPC_WIRE_SIZE * 2; i++) {
        if (hex[i] != ' ' && hex[i] != '\n') clean[ci++] = hex[i];
    }
    clean[ci] = '\0';
    if (ci != UTPC_WIRE_SIZE * 2) {
        fprintf(stderr, "Expected %d hex chars, got %d\n",
                UTPC_WIRE_SIZE * 2, ci);
        return -1;
    }
    uint8_t wire[UTPC_WIRE_SIZE];
    for (int i = 0; i < UTPC_WIRE_SIZE; i++) {
        unsigned int byte;
        if (sscanf(clean + 2*i, "%02x", &byte) != 1) {
            fprintf(stderr, "Bad hex at position %d\n", i);
            return -1;
        }
        wire[i] = (uint8_t)byte;
    }
    utpc_decode(wire, cert);
    return 0;
}

static int verify_file(const char *path) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); return 1; }

    uint8_t wire[UTPC_WIRE_SIZE];
    uint64_t n_total = 0, n_valid = 0, n_invalid = 0;

    while (fread(wire, 1, UTPC_WIRE_SIZE, f) == UTPC_WIRE_SIZE) {
        utpc_t cert;
        utpc_decode(wire, &cert);
        utpc_verify_result_t r = utpc_verify_fast(&cert);
        n_total++;
        if (r.valid) n_valid++;
        else {
            n_invalid++;
            printf("INVALID: ");
            print_result(&cert, &r, 1);
        }
    }
    fclose(f);

    printf("Verified %lu certificates: %lu valid, %lu invalid\n",
           n_total, n_valid, n_invalid);
    return n_invalid > 0 ? 1 : 0;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage:\n");
        printf("  %s <n> <phi> <mu> <J> <dr> <sg>\n", argv[0]);
        printf("  %s --hex <27-byte-hex>\n", argv[0]);
        printf("  %s --file <certs.bin>\n", argv[0]);
        return 1;
    }

    utpc_t cert;

    if (strcmp(argv[1], "--hex") == 0) {
        if (argc < 3) { fprintf(stderr, "--hex requires an argument\n"); return 1; }
        if (parse_hex(argv[2], &cert) != 0) return 1;
    } else if (strcmp(argv[1], "--file") == 0) {
        if (argc < 3) { fprintf(stderr, "--file requires an argument\n"); return 1; }
        return verify_file(argv[2]);
    } else {
        if (argc < 7) {
            fprintf(stderr, "Need 6 values: n phi mu J dr sg\n");
            return 1;
        }
        cert.n   = strtoull(argv[1], NULL, 10);
        cert.phi = strtoull(argv[2], NULL, 10);
        cert.mu  = (uint8_t)atoi(argv[3]);
        cert.J   = strtoull(argv[4], NULL, 10);
        cert.dr  = (uint8_t)atoi(argv[5]);
        cert.sg  = (uint8_t)atoi(argv[6]);
    }

    utpc_verify_result_t r = utpc_verify_fast(&cert);
    print_result(&cert, &r, 1);
    return r.valid ? 0 : 1;
}
