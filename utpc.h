/*
 * utpc.h — Untouchable Totient Prime Certificate
 *
 * Wire format (27 bytes, big-endian):
 *   n    uint64  8 bytes  untouchable number
 *   phi  uint64  8 bytes  euler totient phi(n)
 *   mu   uint8   1 byte   parity divisor (2 or 4)
 *   J    uint64  8 bytes  derived prime J = phi/mu + 1
 *   dr   uint8   1 byte   digital root F10(phi)
 *   sg   uint8   1 byte   sophie germain flag (1 if 2J-1 prime)
 *
 * Construction:
 *   phi = euler_totient(n)
 *   dr  = digital_root(phi)
 *   mu  = (dr >= 5) ? 2 : 4
 *   J   = phi/mu + 1           [must be prime]
 *   sg  = is_prime(2*J - 1)
 *
 * Permanent exclusion (Theorem 1):
 *   p ≡ 7 (mod 9)  =>  no valid UTPC exists for p > 7
 */

#ifndef UTPC_H
#define UTPC_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define UTPC_WIRE_SIZE  27
#define UTPC_VERSION    1

/* ─── Core certificate struct ─────────────────────────────────────────── */

typedef struct {
    uint64_t n;     /* untouchable number                        */
    uint64_t phi;   /* euler totient phi(n)                      */
    uint8_t  mu;    /* parity divisor: 2 or 4                   */
    uint64_t J;     /* derived prime J = phi/mu + 1             */
    uint8_t  dr;    /* digital root of phi                       */
    uint8_t  sg;    /* 1 if 2J-1 is also prime (SG pair)        */
} utpc_t;

/* ─── Verification result ─────────────────────────────────────────────── */

typedef struct {
    int phi_correct;      /* step 1+2: phi == totient(n)        */
    int dr_correct;       /* step 3:   dr == digital_root(phi)  */
    int mu_correct;       /* step 4:   mu consistent with dr    */
    int divisible;        /* step 5:   phi % mu == 0            */
    int J_correct;        /* step 6:   J == phi/mu + 1          */
    int J_prime;          /* step 7:   J is prime               */
    int n_untouchable;    /* step 8:   n is untouchable (slow)  */
    int sg_correct;       /* step 9:   sg flag verified         */
    int valid;            /* all checks passed                  */
} utpc_verify_result_t;

/* ─── Arithmetic helpers ──────────────────────────────────────────────── */

static inline int utpc_digital_root(uint64_t x) {
    if (x == 0) return 0;
    uint64_t r = x % 9;
    return (int)(r == 0 ? 9 : r);
}

static inline uint64_t utpc_modmul(uint64_t a, uint64_t b, uint64_t mod) {
    return (uint64_t)((__uint128_t)a * b % mod);
}

static inline uint64_t utpc_modpow(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = utpc_modmul(result, base, mod);
        base = utpc_modmul(base, base, mod);
        exp >>= 1;
    }
    return result;
}

/* Deterministic Miller-Rabin: correct for all n < 3,317,044,064,679,887,385,961,981 */
static inline int utpc_is_prime(uint64_t n) {
    if (n < 2)  return 0;
    if (n == 2 || n == 3 || n == 5 || n == 7 || n == 11 ||
        n == 13 || n == 17 || n == 19 || n == 23) return 1;
    if (n % 2 == 0 || n % 3 == 0 || n % 5 == 0) return 0;

    uint64_t d = n - 1;
    int s = 0;
    while ((d & 1) == 0) { d >>= 1; s++; }

    static const uint64_t witnesses[] = {
        2, 325, 9375, 28178, 450775, 9780504, 1795265022ULL
    };
    for (int i = 0; i < 7; i++) {
        uint64_t a = witnesses[i] % n;
        if (a == 0) continue;
        uint64_t x = utpc_modpow(a, d, n);
        if (x == 1 || x == n - 1) continue;
        int composite = 1;
        for (int r = 1; r < s; r++) {
            x = utpc_modmul(x, x, n);
            if (x == n - 1) { composite = 0; break; }
        }
        if (composite) return 0;
    }
    return 1;
}

/* Euler totient via trial division */
static inline uint64_t utpc_totient(uint64_t n) {
    uint64_t result = n;
    uint64_t m = n;
    for (uint64_t p = 2; p * p <= m; p++) {
        if (m % p == 0) {
            while (m % p == 0) m /= p;
            result -= result / p;
        }
    }
    if (m > 1) result -= result / m;
    return result;
}

/*
 * Mod-9 fast reject (Theorem 1):
 * Returns 1 if p CANNOT be a J-value (permanent exclusion).
 * p ≡ 7 (mod 9) => blocked (dr(2(p-1))=3, dr(4(p-1))=6, both wrong direction).
 * Call before computing totient to skip ~16.7% of candidates immediately.
 */
static inline int utpc_mod9_blocked(uint64_t p) {
    return (p % 9 == 7);
}

/*
 * Determine mu from digital root.
 * Returns 2 if dr >= 5, 4 if dr < 5, 0 if phi % mu != 0 (invalid).
 */
static inline uint8_t utpc_mu_from_dr(int dr, uint64_t phi) {
    uint8_t mu = (dr >= 5) ? 2 : 4;
    return (phi % mu == 0) ? mu : 0;
}

/* ─── Serialization ───────────────────────────────────────────────────── */

/* Serialize UTPC to 27-byte big-endian wire format */
static inline void utpc_encode(const utpc_t *c, uint8_t out[UTPC_WIRE_SIZE]) {
    uint64_t n   = c->n;
    uint64_t phi = c->phi;
    uint64_t J   = c->J;
    for (int i = 7; i >= 0; i--) { out[i]    = (uint8_t)(n   & 0xFF); n   >>= 8; }
    for (int i = 7; i >= 0; i--) { out[8+i]  = (uint8_t)(phi & 0xFF); phi >>= 8; }
    out[16] = c->mu;
    for (int i = 7; i >= 0; i--) { out[17+i] = (uint8_t)(J   & 0xFF); J   >>= 8; }
    out[25] = c->dr;
    out[26] = c->sg;
}

/* Deserialize 27-byte big-endian wire format to UTPC */
static inline void utpc_decode(const uint8_t in[UTPC_WIRE_SIZE], utpc_t *c) {
    c->n = 0; for (int i=0;i<8;i++) c->n   = (c->n   << 8) | in[i];
    c->phi=0; for (int i=0;i<8;i++) c->phi = (c->phi  << 8) | in[8+i];
    c->mu = in[16];
    c->J = 0; for (int i=0;i<8;i++) c->J   = (c->J   << 8) | in[17+i];
    c->dr = in[25];
    c->sg = in[26];
}

/* Print certificate as hex string */
static inline void utpc_print_hex(const utpc_t *c) {
    uint8_t buf[UTPC_WIRE_SIZE];
    utpc_encode(c, buf);
    for (int i = 0; i < UTPC_WIRE_SIZE; i++) printf("%02x", buf[i]);
}

/* Print certificate in human-readable form */
static inline void utpc_print(const utpc_t *c) {
    printf("UTPC { n=%lu, phi=%lu, mu=%u, J=%lu, dr=%u, sg=%u }",
           c->n, c->phi, c->mu, c->J, c->dr, c->sg);
}

/* ─── Verification (steps 1-7, 9 — fast; step 8 requires sieve) ──────── */

static inline utpc_verify_result_t utpc_verify_fast(const utpc_t *c) {
    utpc_verify_result_t r = {0};

    /* Step 1+2: recompute phi and compare */
    uint64_t phi_check = utpc_totient(c->n);
    r.phi_correct = (phi_check == c->phi);

    /* Step 3: digital root */
    int dr_check = utpc_digital_root(c->phi);
    r.dr_correct = (dr_check == c->dr);

    /* Step 4: mu from dr */
    uint8_t mu_expected = (c->dr >= 5) ? 2 : 4;
    r.mu_correct = (c->mu == mu_expected);

    /* Step 5: divisibility */
    r.divisible = (c->phi % c->mu == 0);

    /* Step 6: J formula */
    uint64_t J_check = c->phi / c->mu + 1;
    r.J_correct = (J_check == c->J);

    /* Step 7: J prime */
    r.J_prime = utpc_is_prime(c->J);

    /* Step 8: skip (requires sieve) */
    r.n_untouchable = -1; /* -1 = not checked */

    /* Step 9: SG flag */
    int sg_actual = utpc_is_prime(2 * c->J - 1);
    r.sg_correct = (sg_actual == c->sg);

    r.valid = r.phi_correct && r.dr_correct && r.mu_correct &&
              r.divisible   && r.J_correct  && r.J_prime    &&
              r.sg_correct;
    return r;
}

#endif /* UTPC_H */
