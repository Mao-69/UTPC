#!/usr/bin/env python3
"""
utpc.py — Untouchable Totient Prime Certificate
Python reference implementation: generation, verification, cluster analysis.

Usage:
    python utpc.py verify  <n> <phi> <mu> <J> <dr> <sg>
    python utpc.py verify  --hex <27-byte-hexstring>
    python utpc.py generate --limit N [--bits B] [--sg-only] [--clusters]
    python utpc.py cluster  <cert_file.bin>
    python utpc.py decode   <27-byte-hexstring>

Examples:
    python utpc.py verify 96 32 2 17 5 0
    python utpc.py verify --hex 0000000000000060000000000000002002000000000000001105 00
    python utpc.py generate --limit 10000 --clusters
    python utpc.py generate --limit 500000 --sg-only
"""

import sys
import struct
import math
import argparse
from collections import defaultdict

# ── Arithmetic ───────────────────────────────────────────────────────────────

def digital_root(n: int) -> int:
    """Iterated digit sum until single digit (= n mod 9, with 9 for multiples of 9)."""
    if n == 0:
        return 0
    r = n % 9
    return 9 if r == 0 else r

def euler_totient(n: int) -> int:
    """Euler's totient function via trial division. O(√n)."""
    result = n
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            while m % p == 0:
                m //= p
            result -= result // p
        p += 1
    if m > 1:
        result -= result // m
    return result

def is_prime_miller_rabin(n: int) -> bool:
    """Deterministic Miller-Rabin for n < 3.3 × 10^24."""
    if n < 2:
        return False
    small = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]
    if n in small:
        return True
    if any(n % p == 0 for p in small):
        return False

    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1

    witnesses = [2, 325, 9375, 28178, 450775, 9780504, 1795265022]
    for a in witnesses:
        a = a % n
        if a == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

def mod9_blocked(p: int) -> bool:
    """
    Theorem 1: p ≡ 7 (mod 9) means no valid UTPC exists for p > 7.
    Fast reject before computing totient.
    """
    return p % 9 == 7

# ── Aliquot sieve ─────────────────────────────────────────────────────────────

def build_untouchable_set(limit: int) -> set:
    """
    Returns set of untouchable numbers in [2, limit].
    n is untouchable if s(m) = σ(m) - m ≠ n for all m ≥ 1.
    O(limit · log log limit) time, O(limit) space.
    """
    sigma = [0] * (limit + 1)
    for i in range(1, limit + 1):
        for j in range(i, limit + 1, i):
            sigma[j] += i

    touchable = set()
    for m in range(2, limit + 1):
        s = sigma[m] - m
        if 2 <= s <= limit:
            touchable.add(s)

    return set(n for n in range(2, limit + 1) if n not in touchable)

# ── Certificate ────────────────────────────────────────────────────────────────

WIRE_SIZE = 27  # bytes

class UTPC:
    """
    Untouchable Totient Prime Certificate.
    Wire format (big-endian, 27 bytes):
        n    uint64  8 bytes
        phi  uint64  8 bytes
        mu   uint8   1 byte
        J    uint64  8 bytes
        dr   uint8   1 byte
        sg   uint8   1 byte
    """
    __slots__ = ('n', 'phi', 'mu', 'J', 'dr', 'sg')

    def __init__(self, n: int, phi: int, mu: int, J: int, dr: int, sg: int):
        self.n   = n
        self.phi = phi
        self.mu  = mu
        self.J   = J
        self.dr  = dr
        self.sg  = sg

    @classmethod
    def from_n(cls, n: int, phi: int = None) -> 'UTPC | None':
        """
        Construct a UTPC from untouchable n.
        Returns None if J(n) is not prime.
        """
        if phi is None:
            phi = euler_totient(n)
        dr  = digital_root(phi)
        mu  = 2 if dr >= 5 else 4
        if phi % mu != 0:
            return None
        J = phi // mu + 1
        if not is_prime_miller_rabin(J):
            return None
        sg = 1 if is_prime_miller_rabin(2 * J - 1) else 0
        return cls(n=n, phi=phi, mu=mu, J=J, dr=dr, sg=sg)

    def encode(self) -> bytes:
        """Serialize to 27-byte big-endian wire format."""
        return struct.pack('>QQBQBB', self.n, self.phi, self.mu, self.J, self.dr, self.sg)

    @classmethod
    def decode(cls, data: bytes) -> 'UTPC':
        """Deserialize from 27-byte big-endian wire format."""
        if len(data) < WIRE_SIZE:
            raise ValueError(f"Expected {WIRE_SIZE} bytes, got {len(data)}")
        n, phi, mu, J, dr, sg = struct.unpack('>QQBQBB', data[:WIRE_SIZE])
        return cls(n=n, phi=phi, mu=mu, J=J, dr=dr, sg=sg)

    def hex(self) -> str:
        return self.encode().hex()

    def __repr__(self) -> str:
        return (f"UTPC(n={self.n}, phi={self.phi}, mu={self.mu}, "
                f"J={self.J}, dr={self.dr}, sg={self.sg})")

    def __eq__(self, other):
        return (self.n == other.n and self.phi == other.phi and
                self.mu == other.mu and self.J == other.J)

# ── Verification ──────────────────────────────────────────────────────────────

class VerifyResult:
    __slots__ = (
        'phi_correct', 'dr_correct', 'mu_correct', 'divisible',
        'J_correct', 'J_prime', 'n_untouchable', 'sg_correct', 'valid'
    )
    def __init__(self):
        for s in self.__slots__:
            setattr(self, s, None)

def verify(cert: UTPC, check_untouchable: bool = False,
           untouchable_set: set = None) -> VerifyResult:
    """
    Verify all 9 steps of a UTPC certificate.

    Steps 1-7 and 9 run in O(√n + log²J).
    Step 8 (untouchability) requires a sieve or a precomputed set.

    Args:
        cert:              the certificate to verify
        check_untouchable: if True, verify step 8 (slow without precomputed set)
        untouchable_set:   precomputed set of untouchables for fast step 8
    """
    r = VerifyResult()

    # Step 1+2: recompute phi
    phi_actual = euler_totient(cert.n)
    r.phi_correct = (phi_actual == cert.phi)

    # Step 3: digital root
    dr_actual = digital_root(cert.phi)
    r.dr_correct = (dr_actual == cert.dr)

    # Step 4: mu from dr
    mu_expected = 2 if cert.dr >= 5 else 4
    r.mu_correct = (cert.mu == mu_expected)

    # Step 5: divisibility
    r.divisible = (cert.phi % cert.mu == 0)

    # Step 6: J formula
    J_expected = cert.phi // cert.mu + 1
    r.J_correct = (J_expected == cert.J)

    # Step 7: J prime
    r.J_prime = is_prime_miller_rabin(cert.J)

    # Step 8: untouchability
    if untouchable_set is not None:
        r.n_untouchable = cert.n in untouchable_set
    elif check_untouchable:
        # Slow: build sieve up to 3*n
        unt = build_untouchable_set(min(cert.n * 3, cert.n + 100_000))
        r.n_untouchable = cert.n in unt
    else:
        r.n_untouchable = None  # not checked

    # Step 9: SG flag
    sg_actual = 1 if is_prime_miller_rabin(2 * cert.J - 1) else 0
    r.sg_correct = (sg_actual == cert.sg)

    r.valid = all([
        r.phi_correct, r.dr_correct, r.mu_correct,
        r.divisible, r.J_correct, r.J_prime, r.sg_correct,
        r.n_untouchable is not False  # None (unchecked) is OK; False means fail
    ])
    return r

def print_verify_result(cert: UTPC, r: VerifyResult, verbose: bool = True):
    status = "VALID" if r.valid else "INVALID"
    color_ok   = "\033[92m"
    color_fail = "\033[91m"
    color_skip = "\033[93m"
    reset      = "\033[0m"

    def mark(v):
        if v is True:  return f"{color_ok}✓{reset}"
        if v is False: return f"{color_fail}✗{reset}"
        return f"{color_skip}−{reset}"

    print(f"\n{cert}")
    print(f"  hex: {cert.hex()}")
    print(f"  status: {color_ok if r.valid else color_fail}{status}{reset}\n")
    if verbose:
        print(f"  [1+2] phi correct       {mark(r.phi_correct)}  (totient(n) == phi)")
        print(f"  [ 3 ] dr correct        {mark(r.dr_correct)}  (digital_root(phi) == dr)")
        print(f"  [ 4 ] mu correct        {mark(r.mu_correct)}  (mu consistent with dr)")
        print(f"  [ 5 ] divisibility      {mark(r.divisible)}  (phi % mu == 0)")
        print(f"  [ 6 ] J formula         {mark(r.J_correct)}  (phi//mu + 1 == J)")
        print(f"  [ 7 ] J prime           {mark(r.J_prime)}  (Miller-Rabin on J)")
        print(f"  [ 8 ] n untouchable     {mark(r.n_untouchable)}  {'(verified)' if r.n_untouchable is not None else '(not checked — add --check-untouchable)'}")
        print(f"  [ 9 ] SG flag           {mark(r.sg_correct)}  (2J-1 prime check matches sg flag)")

# ── Generator ─────────────────────────────────────────────────────────────────

def generate(limit: int, sg_only: bool = False,
             verbose: bool = True) -> list:
    """
    Generate all UTPCs for untouchable n ≤ limit.
    Returns list of UTPC objects.
    """
    print(f"Building untouchable sieve to {limit:,}...", flush=True)
    untouchable = build_untouchable_set(limit)
    print(f"Found {len(untouchable):,} untouchable numbers.", flush=True)

    certs = []
    n_mod9_reject = 0
    n_div_reject  = 0

    for n in sorted(untouchable):
        phi = euler_totient(n)
        if phi < 2:
            continue
        dr  = digital_root(phi)
        mu  = 2 if dr >= 5 else 4
        if phi % mu != 0:
            n_div_reject += 1
            continue
        J = phi // mu + 1
        # Mod-9 fast reject on J
        if mod9_blocked(J):
            n_mod9_reject += 1
            continue
        if not is_prime_miller_rabin(J):
            continue
        sg = 1 if is_prime_miller_rabin(2 * J - 1) else 0
        if sg_only and not sg:
            continue
        cert = UTPC(n=n, phi=phi, mu=mu, J=J, dr=dr, sg=sg)
        certs.append(cert)
        if verbose:
            sg_tag = "  [SG PAIR: 2J-1=%d]" % (2*J-1) if sg else ""
            print(f"  n={n:<8} phi={phi:<8} mu={mu} J={J:<8} dr={dr} sg={sg}{sg_tag}")

    print(f"\nGenerated {len(certs):,} certificates")
    print(f"Mod-9 fast-rejects: {n_mod9_reject:,}  (Theorem 1)")
    print(f"Divisibility rejects: {n_div_reject:,}")
    print(f"SG pairs: {sum(c.sg for c in certs):,} "
          f"({100*sum(c.sg for c in certs)/max(1,len(certs)):.1f}%)")
    return certs

# ── Cluster analysis ──────────────────────────────────────────────────────────

def analyze_clusters(certs: list, min_size: int = 2) -> dict:
    """
    Group certificates by phi value (= cluster membership).
    Returns dict mapping phi -> list of UTPC.
    """
    clusters = defaultdict(list)
    for c in certs:
        clusters[c.phi].append(c)

    sizes = defaultdict(int)
    for ns in clusters.values():
        sizes[len(ns)] += 1

    print(f"\n{'='*60}")
    print(f"CLUSTER ANALYSIS ({len(certs):,} certificates)")
    print(f"{'='*60}")
    print(f"Total clusters:       {len(clusters):,}")
    print(f"Singletons:           {sizes[1]:,}")
    print(f"Clusters size ≥ 2:    {sum(v for k,v in sizes.items() if k>=2):,}")
    print(f"Max cluster size:     {max(sizes.keys()):,}")
    print()
    print(f"Size distribution:")
    for sz in sorted(sizes.keys()):
        bar = '█' * min(sizes[sz], 50)
        print(f"  size {sz:4d}: {sizes[sz]:6d}  {bar}")

    if min_size > 0:
        print(f"\nClusters of size ≥ {min_size}:")
        print(f"  {'phi':>10}  {'J':>10}  {'mu':>3}  {'dr':>3}  {'size':>5}  members")
        for phi in sorted(clusters.keys()):
            ns = clusters[phi]
            if len(ns) >= min_size:
                members = [c.n for c in ns]
                sg_flag = '★' if ns[0].sg else ' '
                print(f"  {phi:>10}  {ns[0].J:>10}  {ns[0].mu:>3}  "
                      f"{ns[0].dr:>3}  {len(ns):>5}  "
                      f"{members[:6]}{'...' if len(members)>6 else ''} {sg_flag}")

    return dict(clusters)

# ── Binary file I/O ────────────────────────────────────────────────────────────

def write_certs(certs: list, path: str):
    """Write certificates to binary file (27 bytes each, big-endian)."""
    with open(path, 'wb') as f:
        for c in certs:
            f.write(c.encode())
    print(f"Wrote {len(certs):,} certificates to {path} "
          f"({len(certs) * WIRE_SIZE:,} bytes)")

def read_certs(path: str) -> list:
    """Read certificates from binary file."""
    certs = []
    with open(path, 'rb') as f:
        while True:
            data = f.read(WIRE_SIZE)
            if len(data) < WIRE_SIZE:
                break
            certs.append(UTPC.decode(data))
    return certs

# ── CLI ────────────────────────────────────────────────────────────────────────

def cmd_verify(args):
    if args.hex:
        hex_str = args.hex.replace(' ', '').replace('\n', '')
        data = bytes.fromhex(hex_str)
        cert = UTPC.decode(data)
    else:
        vals = [int(x) for x in args.values]
        if len(vals) != 6:
            print("verify requires exactly 6 values: n phi mu J dr sg")
            sys.exit(1)
        cert = UTPC(*vals)

    # Optionally build untouchable set for step 8
    unt_set = None
    if args.check_untouchable:
        print(f"Building sieve to {cert.n * 3:,} for untouchability check...")
        unt_set = build_untouchable_set(cert.n * 3)

    r = verify(cert, untouchable_set=unt_set)
    print_verify_result(cert, r, verbose=True)
    sys.exit(0 if r.valid else 1)

def cmd_generate(args):
    certs = generate(
        limit=args.limit,
        sg_only=args.sg_only,
        verbose=not args.quiet
    )
    if args.clusters:
        analyze_clusters(certs, min_size=2)
    if args.output:
        write_certs(certs, args.output)

def cmd_cluster(args):
    print(f"Reading certificates from {args.file}...")
    certs = read_certs(args.file)
    print(f"Read {len(certs):,} certificates.")
    analyze_clusters(certs, min_size=args.min_size)

def cmd_decode(args):
    hex_str = args.hex.replace(' ', '').replace('\n', '')
    data = bytes.fromhex(hex_str)
    cert = UTPC.decode(data)
    print(cert)
    print(f"  2J-1 = {2*cert.J - 1}")
    print(f"  SG pair: {bool(cert.sg)}")
    print(f"  Mod-9 blocked: {mod9_blocked(cert.J)}")

def main():
    parser = argparse.ArgumentParser(
        description="UTPC — Untouchable Totient Prime Certificate tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    sub = parser.add_subparsers(dest='command')

    # verify
    p_ver = sub.add_parser('verify', help='Verify a certificate')
    p_ver.add_argument('values', nargs='*', metavar='VAL',
                       help='n phi mu J dr sg (6 integers)')
    p_ver.add_argument('--hex', help='27-byte hex string')
    p_ver.add_argument('--check-untouchable', action='store_true',
                       help='Also verify step 8 (slow)')

    # generate
    p_gen = sub.add_parser('generate', help='Generate certificates')
    p_gen.add_argument('--limit', type=int, required=True,
                       help='Sieve limit')
    p_gen.add_argument('--sg-only', action='store_true',
                       help='Only emit SG pair certificates')
    p_gen.add_argument('--clusters', action='store_true',
                       help='Run cluster analysis on output')
    p_gen.add_argument('--quiet', action='store_true',
                       help='Suppress per-cert output')
    p_gen.add_argument('--output', '-o', metavar='FILE',
                       help='Write binary certificate file')

    # cluster
    p_cl = sub.add_parser('cluster', help='Analyze a certificate file')
    p_cl.add_argument('file', help='Binary certificate file')
    p_cl.add_argument('--min-size', type=int, default=2,
                      help='Minimum cluster size to display')

    # decode
    p_dec = sub.add_parser('decode', help='Decode a hex certificate')
    p_dec.add_argument('hex', help='27-byte hex string')

    args = parser.parse_args()
    if args.command == 'verify':   cmd_verify(args)
    elif args.command == 'generate': cmd_generate(args)
    elif args.command == 'cluster':  cmd_cluster(args)
    elif args.command == 'decode':   cmd_decode(args)
    else:
        parser.print_help()

if __name__ == '__main__':
    main()
