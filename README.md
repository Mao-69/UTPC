# UTPC — Untouchable Totient Prime Certificate

A prime generation system that produces Sophie Germain primes at **2.41× the
natural rate**, each with a compact 27-byte certificate binding the prime to
its arithmetic lineage through an untouchable number.

---

## Construction

Let `n` be an **untouchable number** (a positive integer not expressible as
the sum of proper divisors of any integer). Define:

```
φ(n)      = Euler totient
dr        = F₁₀(φ(n))  = digital root (iterated digit sum)
μ         = 2 if dr ≥ 5, else 4
J(n)      = φ(n)/μ + 1
```

When `J(n)` is prime, `(n, φ, μ, J)` is a valid **UTPC**.  
When `2J − 1` is also prime, the certificate carries an **SG (Sophie Germain) flag**.

**Theorem 1 (proven):** Any prime `p ≡ 7 (mod 9)` cannot be a J-value.  
This permanently excludes exactly 1/6 of all primes.

---

## Certificate format (27 bytes, big-endian)

| Field | Type   | Bytes | Description                       |
|-------|--------|-------|-----------------------------------|
| `n`   | uint64 | 8     | Untouchable number                |
| `phi` | uint64 | 8     | Euler totient φ(n)                |
| `mu`  | uint8  | 1     | Parity divisor: 2 or 4            |
| `J`   | uint64 | 8     | Derived prime J = φ/μ + 1         |
| `dr`  | uint8  | 1     | Digital root of φ                 |
| `sg`  | uint8  | 1     | 1 if 2J−1 is also prime           |

---

## Build

```bash
make all          # builds utpc_sieve and utpc_verify_cli
make test         # runs built-in tests
```

Requirements: `gcc`, `pthread`, Python 3.8+ (for `utpc.py`)

---

## Usage

### C sieve (fast, multi-threaded)

```bash
# Sieve to 1 million, 8 threads
./utpc_sieve 1000000 8

# SG pairs only, quiet mode
./utpc_sieve --sg-only --quiet 500000000 $(nproc)

# With cluster analysis
./utpc_sieve --clusters 10000000 4

# Custom output file
./utpc_sieve --output my_certs.bin --log my_log.txt 1000000 8
```

Output files:
- `utpc_certs.bin` — binary certificates (27 bytes each)
- `utpc_results.txt` — human-readable log
- `utpc_clusters.txt` — cluster table (with `--clusters`)

### C verifier

```bash
# Verify by fields
./utpc_verify_cli 96 32 2 17 5 0

# Verify by hex
./utpc_verify_cli --hex 0000000000000060000000000000002002000000000000001105 00

# Verify all certs in a file
./utpc_verify_cli --file utpc_certs.bin
```

### Python reference implementation

```bash
# Verify a certificate
python3 utpc.py verify 96 32 2 17 5 0
python3 utpc.py verify --check-untouchable 96 32 2 17 5 0

# Generate up to limit, with cluster analysis
python3 utpc.py generate --limit 10000 --clusters
python3 utpc.py generate --limit 500000 --sg-only --output sg_certs.bin

# Analyze an existing certificate file
python3 utpc.py cluster utpc_certs.bin --min-size 3

# Decode a hex certificate
python3 utpc.py decode 0000000000000060000000000000002002000000000000001105 00
```
```
python3 utpc.py verify --check-untouchable 146 72 2 37 9 1
Building sieve to 438 for untouchability check...

UTPC(n=146, phi=72, mu=2, J=37, dr=9, sg=1)
  hex: 000000000000009200000000000000480200000000000000250901
  status: VALID

  [1+2] phi correct       ✓  (totient(n) == phi)
  [ 3 ] dr correct        ✓  (digital_root(phi) == dr)
  [ 4 ] mu correct        ✓  (mu consistent with dr)
  [ 5 ] divisibility      ✓  (phi % mu == 0)
  [ 6 ] J formula         ✓  (phi//mu + 1 == J)
  [ 7 ] J prime           ✓  (Miller-Rabin on J)
  [ 8 ] n untouchable     ✓  (verified)
  [ 9 ] SG flag           ✓  (2J-1 prime check matches sg flag)
```
```
python3 ./utpc.py decode 000000000000009200000000000000480200000000000000250901
UTPC(n=146, phi=72, mu=2, J=37, dr=9, sg=1)
  2J-1 = 73
  SG pair: True
  Mod-9 blocked: False
```
### Extract first certificate as hex (for inspection)
```
python3 -c '
with open("utpc_certs.bin", "rb") as f:
    data = f.read(27)
    print(data.hex())
'
000000000000000500000000000000040400000000000000020401
```
### Use the Python API directly (cleanest for single generation):
```
from utpc import UTPC

# Generate one certificate from an untouchable number
cert = UTPC.from_n(96)
print(cert)
print("Hex (27 bytes):", cert.hex())

# Another example
cert2 = UTPC.from_n(146)
print(cert2)
print("Hex:", cert2.hex())
```
```
UTPC(n=96, phi=32, mu=2, J=17, dr=5, sg=0)
Hex (27 bytes): 000000000000006000000000000000200200000000000000110500
UTPC(n=146, phi=72, mu=2, J=37, dr=9, sg=1)
Hex: 000000000000009200000000000000480200000000000000250901
```

### Python API

```python
from utpc import UTPC, verify, generate, analyze_clusters, build_untouchable_set

# Build a certificate from n
cert = UTPC.from_n(96)
print(cert)          # UTPC(n=96, phi=32, mu=2, J=17, dr=5, sg=0)
print(cert.hex())    # 27-byte hex

# Verify
unt = build_untouchable_set(300)   # precompute for fast step 8
result = verify(cert, untouchable_set=unt)
print(result.valid)  # True

# Encode / decode
wire = cert.encode()               # 27 bytes
cert2 = UTPC.decode(wire)          # round-trips exactly

# Generate and cluster
certs = generate(limit=10000, sg_only=False)
clusters = analyze_clusters(certs, min_size=2)
```

---

## Sample certificates

| n   | φ(n) | μ | J  | dr | SG | hex (27 bytes)                                         |
|-----|------|---|----|----|----|--------------------------------------------------------|
| 96  | 32   | 2 | 17 | 5  | No | `000000000000006000000000000000200200000000000000110500` |
| 120 | 32   | 2 | 17 | 5  | No | `000000000000007800000000000000200200000000000000110500` |
| 146 | 72   | 2 | 37 | 9  | Yes| `000000000000009200000000000000480200000000000000250901` |
| 216 | 72   | 2 | 37 | 9  | Yes| `00000000000000d800000000000000480200000000000000250901` |
| 292 | 144  | 2 | 73 | 9  | No | `000000000000012400000000000000900200000000000000490900` |
| 304 | 144  | 2 | 73 | 9  | No | `000000000000013000000000000000900200000000000000490900` |

---

## Key results (sieve to 500,000,000)

| Metric                     | Value                    |
|----------------------------|--------------------------|
| Qualifying J-prime pairs   | 57,498,184               |
| μ=2 / μ=4 ratio            | 3.148  (theory: 4.0)     |
| SG pairs                   | 7,171,025  (12.47%)      |
| SG enrichment over baseline| **2.41×**                |
| Largest J-prime            | 249,999,661              |
| Max cluster size           | 306 (φ=207360, J=103681) |

---

## Cryptographic applications

| Use case            | Role of J                    | Standard           |
|---------------------|------------------------------|--------------------|
| Diffie–Hellman      | p = 2J−1 safe prime          | RFC 7919           |
| RSA strong primes   | p−1 has factor J ≈ p/2       | FIPS 186-4 §B.3.3  |
| DSA parameters      | q = J, p = 2J−1              | FIPS 186-4 §A.1    |
| Threshold / Paillier| Cluster = distributed key gen| Paillier 1999      |
| HSM attestation     | 27-byte cert for constrained devices | FIPS 140-3  |
| Post-quantum hybrid | Classical DH safe prime      | NIST SP 800-56C    |

---

## Related sequences (OEIS)

- A005114 — Untouchable numbers
- A000010 — Euler's totient function  
- A005384 — Sophie Germain primes
- J-value sequence and cluster size sequence — submitted / pending

---

## Open problems | future work

1. Derive the SG enrichment constant (~1.83× asymptotic) analytically
2. Prove cluster sizes are unbounded
3. Prove the asymptotic J-coverage conjecture (5/6 of all primes)
4. Formal security proof for the cluster DKG protocol
5. Sub-sieve totient inversion algorithm (or proof of hardness)
