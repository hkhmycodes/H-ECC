// exp4_openmp.cpp  (rev 2)
//
// EXP-4: CPU parallel scaling benchmark for the pK representation versus
// Jacobian coordinates on the Weierstrass model E_W, using OpenMP and GMP.
//
// Changes relative to rev 1:
//   1. jac_equal() now compares Jacobian points correctly
//      (X1*Z2^2 = X2*Z1^2, Y1*Z2^3 = Y2*Z1^3). The rev-1 test used the
//      homogeneous rule X1*Z2 = X2*Z1, which is wrong for Jacobian
//      coordinates and made the pre-timing validation fail.
//   2. Validation is no longer tautological: it cross-checks
//      Phi_D([k]P^pK) against [k]Phi_D(P^pK) and checks [q]P^W = O.
//   3. After timing (outside the timed region), EVERY pK output is mapped
//      through Phi_D and compared with the Jacobian output.
//   4. Cross-thread identity is checked on the full canonical outputs
//      (exact equality), not only on a 64-bit checksum of low limbs.
//      mpz_get_ui() ignores the sign, and GMP's % can return negative
//      residues, so values are reduced to [0,P) before hashing.
//   5. The actual OpenMP team size is verified for every T.
//   6. All 5 timing samples are written to the JSON, not only the median.
//
// Build (MSYS2 MINGW64 or Linux):
//   g++ -std=c++17 -O3 -fopenmp -o exp4_openmp exp4_openmp.cpp -lgmpxx -lgmp
//
// Run:
//   ./exp4_openmp <N> <seed> [T1 T2 T3 ...]
//
// stdout: one JSON object.  stderr: progress.
// Exit codes: 0 ok, 1 usage, 2 pre-timing validation failed,
//             3 outputs differ across thread counts,
//             4 pK/Jacobian oracle mismatch, 5 OpenMP team size != T.

#include <gmpxx.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>

// ----------------------------------------------------------------
// Constants (test vector of Section 6.5)
// ----------------------------------------------------------------
static const char* P_STR =
    "106839527430202782610735740077697524141755216002708025221579380756122101340723";
static const char* Q_STR =
    "26709881857550695652683935019424381035438804000677006305394845189030525335181";
static const char* G_B_STR =
    "43885198696659819331868805857976446729580581201517994933321256281291666494433";
static const char* G_C_STR =
    "52096118237682591038514212999531378342298973921021687132305677862559506645034";

static mpz_class P, Q_ORDER;
static mpz_class A_W;                      // coefficient of x in E_W
static mpz_class INV2, INV3, INV6;

static const int REPEATS = 5;

struct PkPoint  { mpz_class a, b, c; };
struct JacPoint { mpz_class X, Y, Z; };

// ----------------------------------------------------------------
// Composite pK arithmetic (Sections 4.1 and 4.2 of the manuscript)
// ----------------------------------------------------------------
static inline void pk_add(const PkPoint& A, const PkPoint& B, PkPoint& R) {
    const mpz_class& a = A.a; const mpz_class& b = A.b; const mpz_class& c = A.c;
    const mpz_class& x = B.a; const mpz_class& y = B.b; const mpz_class& z = B.c;
    mpz_class Fa = (5 * (b*b - c*c) + 2 * a * (7*c - b)) % P;
    mpz_class Fb = (   (c*c - a*a) + 2 * b * (5*a - 7*c)) % P;
    mpz_class Fc = (7 * (a*a - b*b) + 2 * c * (b - 5*a)) % P;
    mpz_class Ga = (5 * (y*y - z*z) + 2 * x * (7*z - y)) % P;
    mpz_class Gb = (   (z*z - x*x) + 2 * y * (5*x - 7*z)) % P;
    mpz_class Gc = (7 * (x*x - y*y) + 2 * z * (y - 5*x)) % P;
    mpz_class A_ = (a*Ga + b*Gb + c*Gc) % P;
    mpz_class B_ = (x*Fa + y*Fb + z*Fc) % P;
    mpz_class R0 = (A_*a - B_*x) % P;
    mpz_class R1 = (A_*b - B_*y) % P;
    mpz_class R2 = (A_*c - B_*z) % P;
    // All reads of A and B are finished; R may alias A.
    R.a = (R1 * R2) % P;
    R.b = (R2 * R0) % P;
    R.c = (R0 * R1) % P;
}

static inline void pk_double(const PkPoint& A, PkPoint& R) {
    const mpz_class& a = A.a; const mpz_class& b = A.b; const mpz_class& c = A.c;
    mpz_class Fa = (5 * (b*b - c*c) + 2 * a * (7*c - b)) % P;
    mpz_class Fb = (   (c*c - a*a) + 2 * b * (5*a - 7*c)) % P;
    mpz_class Fc = (7 * (a*a - b*b) + 2 * c * (b - 5*a)) % P;
    mpz_class B_ = Fb, G_ = Fc;
    mpz_class Ft = (G_*B_*B_ + 7*B_*G_*G_) % P;
    mpz_class H  = ((5*a - 7*c)*G_*G_
                    - 2*(c - 7*b)*B_*G_
                    + (b - 5*a)*B_*B_) % P;
    mpz_class R0 = (-Ft*a) % P;
    mpz_class R1 = (-Ft*b + H*G_) % P;
    mpz_class R2 = (-Ft*c - H*B_) % P;
    R.a = (R1 * R2) % P;
    R.b = (R2 * R0) % P;
    R.c = (R0 * R1) % P;
}

// ----------------------------------------------------------------
// Jacobian arithmetic on E_W (x = X/Z^2, y = Y/Z^3)
// ----------------------------------------------------------------
static inline void jac_add(const JacPoint& A, const JacPoint& B, JacPoint& R) {
    const mpz_class& X1 = A.X; const mpz_class& Y1 = A.Y; const mpz_class& Z1 = A.Z;
    const mpz_class& X2 = B.X; const mpz_class& Y2 = B.Y; const mpz_class& Z2 = B.Z;
    mpz_class Z1Z1 = (Z1*Z1) % P, Z2Z2 = (Z2*Z2) % P;
    mpz_class U1 = (X1*Z2Z2) % P, U2 = (X2*Z1Z1) % P;
    mpz_class S1 = (Y1*Z2 % P * Z2Z2) % P;
    mpz_class S2 = (Y2*Z1 % P * Z1Z1) % P;
    mpz_class H  = (U2 - U1) % P;
    mpz_class I  = (4 * H * H) % P;
    mpz_class J  = (H * I) % P;
    mpz_class r  = (2 * (S2 - S1)) % P;
    mpz_class VV = (U1 * I) % P;
    mpz_class X3 = (r*r - J - 2*VV) % P;
    mpz_class Y3 = (r * (VV - X3) - 2 * S1 * J) % P;
    mpz_class t  = (Z1 + Z2); t = (t*t - Z1Z1 - Z2Z2) % P;
    mpz_class Z3 = (t * H) % P;
    R.X = X3; R.Y = Y3; R.Z = Z3;
}

static inline void jac_dbl(const JacPoint& A, JacPoint& R) {
    const mpz_class& X1 = A.X; const mpz_class& Y1 = A.Y; const mpz_class& Z1 = A.Z;
    mpz_class XX = (X1*X1) % P, YY = (Y1*Y1) % P;
    mpz_class YYYY = (YY*YY) % P, ZZ = (Z1*Z1) % P;
    mpz_class t = (X1 + YY); t = (t*t - XX - YYYY) % P;
    mpz_class S = (2 * t) % P;
    mpz_class M = (3*XX + A_W * ((ZZ*ZZ) % P)) % P;
    mpz_class T = (M*M - 2*S) % P;
    mpz_class t2 = (Y1 + Z1); t2 = (t2*t2 - YY - ZZ) % P;
    mpz_class Y3 = (M * (S - T) - 8 * YYYY) % P;
    R.X = T; R.Y = Y3; R.Z = t2;
}

static inline bool is_zero_mod(const mpz_class& v) {
    mpz_class r = v % P;
    return r == 0;
}

// Correct projective equality for Jacobian coordinates.
static bool jac_equal(const JacPoint& A, const JacPoint& B) {
    bool ia = is_zero_mod(A.Z), ib = is_zero_mod(B.Z);
    if (ia || ib) return ia && ib;
    mpz_class Z1Z1 = (A.Z * A.Z) % P, Z2Z2 = (B.Z * B.Z) % P;
    if (!is_zero_mod(A.X * Z2Z2 - B.X * Z1Z1)) return false;
    mpz_class Z1c = (Z1Z1 * A.Z) % P, Z2c = (Z2Z2 * B.Z) % P;
    if (!is_zero_mod(A.Y * Z2c - B.Y * Z1c)) return false;
    return true;
}

// ----------------------------------------------------------------
// Scalar multiplication (right-to-left double-and-add)
// ----------------------------------------------------------------
static void pk_scalar_mul(const mpz_class& k, const PkPoint& Pt, PkPoint& R) {
    PkPoint D; D.a = 5; D.b = 1; D.c = 7;
    if (k == 0) { R = D; return; }
    PkPoint Q = Pt, acc = D;
    mpz_class kk = k;
    while (kk > 0) {
        if (mpz_odd_p(kk.get_mpz_t())) pk_add(acc, Q, acc);
        pk_double(Q, Q);
        kk >>= 1;
    }
    R = acc;
}

static void jac_scalar_mul(const mpz_class& k, const JacPoint& Pt, JacPoint& R) {
    JacPoint inf; inf.X = 0; inf.Y = 1; inf.Z = 0;
    if (k == 0) { R = inf; return; }
    JacPoint Q = Pt;
    bool have = false;
    JacPoint acc;
    mpz_class kk = k;
    while (kk > 0) {
        if (mpz_odd_p(kk.get_mpz_t())) {
            if (!have) { acc = Q; have = true; }
            else       { jac_add(acc, Q, acc); }
        }
        jac_dbl(Q, Q);
        kk >>= 1;
    }
    R = have ? acc : inf;
}

// ----------------------------------------------------------------
// Affine Weierstrass helper for the birational map
// ----------------------------------------------------------------
static void w_add_affine(const mpz_class& x1, const mpz_class& y1,
                         const mpz_class& x2, const mpz_class& y2,
                         mpz_class& x3, mpz_class& y3, bool& inf) {
    if ((x1 - x2) % P == 0 && (y1 + y2) % P == 0) { inf = true; return; }
    mpz_class lam;
    if (x1 == x2 && y1 == y2) {
        mpz_class num = (3*x1*x1 + A_W) % P;
        mpz_class den = (2*y1) % P;
        mpz_class dinv; mpz_invert(dinv.get_mpz_t(), den.get_mpz_t(), P.get_mpz_t());
        lam = (num * dinv) % P;
    } else {
        mpz_class num = (y2 - y1) % P; if (num < 0) num += P;
        mpz_class den = (x2 - x1) % P; if (den < 0) den += P;
        mpz_class dinv; mpz_invert(dinv.get_mpz_t(), den.get_mpz_t(), P.get_mpz_t());
        lam = (num * dinv) % P;
    }
    x3 = ((lam*lam - x1 - x2) % P + P) % P;
    y3 = ((lam*(x1 - x3) - y1) % P + P) % P;
    inf = false;
}

// ----------------------------------------------------------------
// Phi_D: pK -> E_W (Jacobian, Z = 1)
//
// For (u:v:w) = (5:1:7):
//   a_map = 2, d = -1/6, lambda = -1, shift = 0
//   r = (x-1)/(x+1), s = (y-1)/(y+1), t = r/s, T = t/2
//   X = -T/6,  Y = -s*T*(T-1)/6
//   W_D = (1/3, -1/6)
//   Phi_D(P) = (X,Y) - W_D
// ----------------------------------------------------------------
static bool phi_D_to_W(const PkPoint& p, JacPoint& out) {
    if (is_zero_mod(p.a)) return false;
    mpz_class ai; if (mpz_invert(ai.get_mpz_t(), p.a.get_mpz_t(), P.get_mpz_t()) == 0)
        return false;
    mpz_class x = (p.b * ai) % P;
    mpz_class y = (p.c * ai) % P;

    mpz_class xp1 = (x + 1) % P, yp1 = (y + 1) % P;
    if (xp1 == 0 || yp1 == 0) return false;
    mpz_class xm1 = (x - 1) % P; if (xm1 < 0) xm1 += P;
    mpz_class ym1 = (y - 1) % P; if (ym1 < 0) ym1 += P;
    mpz_class xp1i, yp1i;
    mpz_invert(xp1i.get_mpz_t(), xp1.get_mpz_t(), P.get_mpz_t());
    mpz_invert(yp1i.get_mpz_t(), yp1.get_mpz_t(), P.get_mpz_t());
    mpz_class r = (xm1 * xp1i) % P;
    mpz_class s = (ym1 * yp1i) % P;
    if (s == 0) return false;
    mpz_class si; mpz_invert(si.get_mpz_t(), s.get_mpz_t(), P.get_mpz_t());
    mpz_class t = (r * si) % P;
    mpz_class T = (t * INV2) % P;
    mpz_class XW = (P - (T * INV6) % P) % P;
    mpz_class Tm1 = (T - 1) % P; if (Tm1 < 0) Tm1 += P;
    mpz_class Y = (P - (((s * T) % P * Tm1) % P * INV6) % P) % P;

    mpz_class WDx = INV3;
    mpz_class WDy = (P - INV6) % P;
    mpz_class nWDy = (P - WDy) % P;
    mpz_class X3, Y3; bool inf;
    w_add_affine(XW, Y, WDx, nWDy, X3, Y3, inf);
    if (inf) return false;
    out.X = X3; out.Y = Y3; out.Z = 1;
    return true;
}

// ----------------------------------------------------------------
// RNG
// ----------------------------------------------------------------
static inline void mpz_from_rng(mpz_class& out, std::mt19937_64& rng) {
    std::array<uint64_t, 4> w{ rng(), rng(), rng(), rng() };
    mpz_import(out.get_mpz_t(), 4, -1, sizeof(uint64_t), 0, 0, w.data());
}

// Timing workload only: scalar in [1, q-1]. (Not a uniform sampler.)
static mpz_class rng_scalar_in_subgroup(std::mt19937_64& rng) {
    mpz_class k;
    mpz_from_rng(k, rng);
    k = (k % (Q_ORDER - 1)) + 1;
    return k;
}

// ----------------------------------------------------------------
// Workload generation (once per (N, seed), before any timing)
// ----------------------------------------------------------------
struct Workload {
    std::vector<mpz_class> scalars;
    std::vector<PkPoint>   pk_pts;
    std::vector<JacPoint>  jac_pts;
    unsigned n_retries = 0;
};

static Workload make_workload(unsigned N, unsigned seed) {
    Workload W;
    W.scalars.reserve(N); W.pk_pts.reserve(N); W.jac_pts.reserve(N);
    std::mt19937_64 rng(((uint64_t)seed << 32) ^ 0x9E3779B97F4A7C15ULL);
    PkPoint G; G.a = 1;
    G.b.set_str(G_B_STR, 10);
    G.c.set_str(G_C_STR, 10);
    while (W.pk_pts.size() < N) {
        mpz_class s = rng_scalar_in_subgroup(rng);
        PkPoint Pi;
        pk_scalar_mul(s, G, Pi);
        JacPoint Wi;
        if (!phi_D_to_W(Pi, Wi)) {
            ++W.n_retries;
            continue;
        }
        W.pk_pts.push_back(Pi);
        W.jac_pts.push_back(Wi);
        W.scalars.push_back(rng_scalar_in_subgroup(rng));
    }
    return W;
}

// ----------------------------------------------------------------
// Pre-timing validation (subset)
//   (a) [q]P^W = O                     (point lies in the order-q subgroup)
//   (b) Phi_D([k]P^pK) == [k]Phi_D(P^pK)   (pK arithmetic vs Weierstrass)
// ----------------------------------------------------------------
static bool validate_subset(const Workload& W, unsigned subset_size) {
    unsigned m = std::min<unsigned>(subset_size, (unsigned)W.scalars.size());
    for (unsigned i = 0; i < m; ++i) {
        JacPoint qP;
        jac_scalar_mul(Q_ORDER, W.jac_pts[i], qP);
        if (!is_zero_mod(qP.Z)) {
            std::fprintf(stderr, "[validate] [q]P_%u^W != O\n", i);
            return false;
        }
        PkPoint out_pk;
        pk_scalar_mul(W.scalars[i], W.pk_pts[i], out_pk);
        JacPoint expected_W;
        if (!phi_D_to_W(out_pk, expected_W)) {
            std::fprintf(stderr, "[validate] Phi_D failed on [k_%u]P_%u^pK\n", i, i);
            return false;
        }
        JacPoint got_W;
        jac_scalar_mul(W.scalars[i], W.jac_pts[i], got_W);
        if (!jac_equal(expected_W, got_W)) {
            std::fprintf(stderr, "[validate] [k_%u]P_%u: pK->W and W disagree\n", i, i);
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------
// Canonical flattening and hashing (outside timed region)
// ----------------------------------------------------------------
static inline mpz_class canon(const mpz_class& v) {
    mpz_class r = v % P;
    if (r < 0) r += P;
    return r;
}

static std::vector<mpz_class> flatten(const std::vector<PkPoint>& v) {
    std::vector<mpz_class> f; f.reserve(3 * v.size());
    for (const auto& p : v) { f.push_back(canon(p.a)); f.push_back(canon(p.b)); f.push_back(canon(p.c)); }
    return f;
}
static std::vector<mpz_class> flatten(const std::vector<JacPoint>& v) {
    std::vector<mpz_class> f; f.reserve(3 * v.size());
    for (const auto& p : v) { f.push_back(canon(p.X)); f.push_back(canon(p.Y)); f.push_back(canon(p.Z)); }
    return f;
}

static uint64_t hash_flat(const std::vector<mpz_class>& f) {
    uint64_t h = 1469598103934665603ULL;               // FNV-1a style over all limbs
    for (const auto& x : f) {
        size_t n = mpz_size(x.get_mpz_t());
        for (size_t i = 0; i < n; ++i) {
            h ^= (uint64_t)mpz_getlimbn(x.get_mpz_t(), i);
            h *= 1099511628211ULL;
            h ^= h >> 29;
        }
        h ^= 0xFFULL; h *= 1099511628211ULL;
    }
    return h;
}

// ----------------------------------------------------------------
// Timing kernels. T=1: plain loop, no OpenMP region.  T>1: omp parallel for.
// ----------------------------------------------------------------
struct TimingResult { double median; std::vector<double> samples; };

static TimingResult finish(std::vector<double> s) {
    std::vector<double> sorted = s;
    std::sort(sorted.begin(), sorted.end());
    return { sorted[sorted.size() / 2], s };
}

template <class Pt, class Mul>
static TimingResult run_seq(const std::vector<mpz_class>& ks,
                            const std::vector<Pt>& pts, Mul mul,
                            int repeats, std::vector<Pt>& outs) {
    const size_t N = ks.size();
    outs.assign(N, Pt());
    for (size_t i = 0; i < N; ++i) mul(ks[i], pts[i], outs[i]);          // warm-up
    std::vector<double> s(repeats);
    for (int r = 0; r < repeats; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        for (size_t i = 0; i < N; ++i) mul(ks[i], pts[i], outs[i]);
        auto t1 = std::chrono::steady_clock::now();
        s[r] = std::chrono::duration<double>(t1 - t0).count();
    }
    return finish(s);
}

template <class Pt, class Mul>
static TimingResult run_par(const std::vector<mpz_class>& ks,
                            const std::vector<Pt>& pts, Mul mul,
                            int T, int repeats, std::vector<Pt>& outs) {
    const long long N = (long long)ks.size();
    outs.assign((size_t)N, Pt());
    #pragma omp parallel for schedule(static) num_threads(T)             // warm-up
    for (long long i = 0; i < N; ++i) mul(ks[(size_t)i], pts[(size_t)i], outs[(size_t)i]);
    std::vector<double> s(repeats);
    for (int r = 0; r < repeats; ++r) {
        auto t0 = std::chrono::steady_clock::now();
        #pragma omp parallel for schedule(static) num_threads(T)
        for (long long i = 0; i < N; ++i) mul(ks[(size_t)i], pts[(size_t)i], outs[(size_t)i]);
        auto t1 = std::chrono::steady_clock::now();
        s[r] = std::chrono::duration<double>(t1 - t0).count();
    }
    return finish(s);
}

static int actual_team_size(int T) {
    int n = 0;
    #pragma omp parallel num_threads(T)
    {
        #pragma omp single
        n = omp_get_num_threads();
    }
    return n;
}

// ----------------------------------------------------------------
// Per-representation driver
// ----------------------------------------------------------------
struct Row {
    std::string kind;
    int T, team;
    double median, throughput, speedup, efficiency;
    uint64_t hash;
    bool identical;
    std::vector<double> samples;
};

template <class Pt, class Mul>
static void run_kind(const char* kind, const Workload& W,
                     const std::vector<Pt>& pts, Mul mul,
                     const std::vector<int>& threads, int repeats,
                     std::vector<Pt>& ref_out, std::vector<Row>& rows) {
    const unsigned N = (unsigned)W.scalars.size();

    TimingResult seq = run_seq(W.scalars, pts, mul, repeats, ref_out);
    std::vector<mpz_class> ref_flat = flatten(ref_out);
    uint64_t h_ref = hash_flat(ref_flat);
    std::fprintf(stderr, "[exp4] %-3s N=%u T=1  median=%.6g s  hash=%016llx\n",
                 kind, N, seq.median, (unsigned long long)h_ref);
    rows.push_back({ kind, 1, 1, seq.median, N / seq.median, 1.0, 1.0,
                     h_ref, true, seq.samples });

    for (int T : threads) {
        if (T <= 1 || (unsigned)T > N) continue;
        int team = actual_team_size(T);
        std::vector<Pt> outs;
        TimingResult par = run_par(W.scalars, pts, mul, T, repeats, outs);
        std::vector<mpz_class> flat = flatten(outs);
        bool same = (flat == ref_flat);                 // exact equality
        double sp = seq.median / par.median;
        std::fprintf(stderr,
            "[exp4] %-3s N=%u T=%d (team=%d)  median=%.6g s  speedup=%.3fx  "
            "eff=%.3f  hash=%016llx %s\n",
            kind, N, T, team, par.median, sp, sp / T,
            (unsigned long long)hash_flat(flat), same ? "IDENTICAL" : "MISMATCH");
        rows.push_back({ kind, T, team, par.median, N / par.median, sp, sp / T,
                         hash_flat(flat), same, par.samples });
    }
}

// ----------------------------------------------------------------
// Post-timing oracle: every pK output vs the Jacobian output
// ----------------------------------------------------------------
static void oracle_full(const std::vector<PkPoint>& pk,
                        const std::vector<JacPoint>& jc,
                        unsigned& checked, unsigned& skipped, unsigned& bad) {
    checked = skipped = bad = 0;
    for (size_t i = 0; i < pk.size(); ++i) {
        JacPoint m;
        if (!phi_D_to_W(pk[i], m)) { ++skipped; continue; }
        ++checked;
        if (!jac_equal(m, jc[i])) ++bad;
    }
}

static const char* envs(const char* k) { const char* v = std::getenv(k); return v ? v : ""; }

static void print_rows(const std::vector<Row>& rows, bool& first) {
    for (const auto& r : rows) {
        std::printf("%s    {\"kind\":\"%s\",\"T\":%d,\"team\":%d,"
                    "\"T_median_s\":%.6g,\"throughput_ops\":%.6g,"
                    "\"speedup\":%.6g,\"efficiency\":%.6g,"
                    "\"hash\":\"%016llx\",\"identical_to_T1\":%s,\"samples_s\":[",
                    first ? "" : ",\n", r.kind.c_str(), r.T, r.team,
                    r.median, r.throughput, r.speedup, r.efficiency,
                    (unsigned long long)r.hash, r.identical ? "true" : "false");
        for (size_t i = 0; i < r.samples.size(); ++i)
            std::printf("%s%.6g", i ? "," : "", r.samples[i]);
        std::printf("]}");
        first = false;
    }
}

// ----------------------------------------------------------------
// Main
// ----------------------------------------------------------------
int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s <N> <seed> [T1 T2 T3 ...]\n", argv[0]);
        return 1;
    }
    unsigned N    = (unsigned)std::atoi(argv[1]);
    unsigned seed = (unsigned)std::atoi(argv[2]);
    if (N == 0) { std::fprintf(stderr, "N must be >= 1\n"); return 1; }

    std::vector<int> threads;
    if (argc == 3) threads = { 1, 2, 4, 8 };
    else for (int i = 3; i < argc; ++i) threads.push_back(std::atoi(argv[i]));

    omp_set_dynamic(0);

    P.set_str(P_STR, 10);
    Q_ORDER.set_str(Q_STR, 10);
    mpz_invert(INV2.get_mpz_t(), mpz_class(2).get_mpz_t(), P.get_mpz_t());
    mpz_invert(INV3.get_mpz_t(), mpz_class(3).get_mpz_t(), P.get_mpz_t());
    mpz_invert(INV6.get_mpz_t(), mpz_class(6).get_mpz_t(), P.get_mpz_t());
    mpz_class INV36; mpz_invert(INV36.get_mpz_t(), mpz_class(36).get_mpz_t(), P.get_mpz_t());
    A_W = (P - INV36) % P;

    // --- Workload generation (NOT timed) ---
    auto tg0 = std::chrono::steady_clock::now();
    Workload W = make_workload(N, seed);
    auto tg1 = std::chrono::steady_clock::now();
    double t_gen = std::chrono::duration<double>(tg1 - tg0).count();
    std::fprintf(stderr, "[exp4] workload generation: %.3f s (%u retries)\n",
                 t_gen, W.n_retries);

    // --- Pre-timing validation ---
    unsigned subset = std::min<unsigned>(20u, N);
    bool valid = validate_subset(W, subset);
    if (!valid) {
        std::fprintf(stderr, "[exp4] validation FAILED; aborting\n");
        return 2;
    }
    std::fprintf(stderr, "[exp4] pre-timing validation: OK on %u/%u points\n", subset, N);

    // --- Timed runs ---
    std::vector<Row> rows_pk, rows_jac;
    std::vector<PkPoint>  ref_pk;
    std::vector<JacPoint> ref_jac;
    run_kind("pk", W, W.pk_pts, [](const mpz_class& k, const PkPoint& p, PkPoint& r) {
        pk_scalar_mul(k, p, r); }, threads, REPEATS, ref_pk, rows_pk);
    run_kind("jac", W, W.jac_pts, [](const mpz_class& k, const JacPoint& p, JacPoint& r) {
        jac_scalar_mul(k, p, r); }, threads, REPEATS, ref_jac, rows_jac);

    // --- Post-timing checks (not timed) ---
    unsigned o_checked, o_skipped, o_bad;
    oracle_full(ref_pk, ref_jac, o_checked, o_skipped, o_bad);
    bool oracle_ok = (o_bad == 0 && o_checked > 0);
    std::fprintf(stderr, "[exp4] oracle (all outputs): checked=%u skipped=%u mismatches=%u  %s\n",
                 o_checked, o_skipped, o_bad, oracle_ok ? "OK" : "FAILED");

    bool all_identical = true, team_ok = true;
    for (const auto* rs : { &rows_pk, &rows_jac })
        for (const auto& r : *rs) {
            if (!r.identical) all_identical = false;
            if (r.team != r.T) team_ok = false;
        }

    // --- JSON ---
    std::printf("{\n");
    std::printf("  \"version\": \"rev2\",\n");
    std::printf("  \"n\": %u,\n  \"seed\": %u,\n  \"repeats\": %d,\n", N, seed, REPEATS);
    std::printf("  \"workload_generation_s\": %.6g,\n", t_gen);
    std::printf("  \"validation_subset\": %u,\n  \"validation_passed\": %s,\n",
                subset, valid ? "true" : "false");
    std::printf("  \"oracle_checked\": %u,\n  \"oracle_skipped\": %u,\n"
                "  \"oracle_mismatches\": %u,\n  \"oracle_passed\": %s,\n",
                o_checked, o_skipped, o_bad, oracle_ok ? "true" : "false");
    std::printf("  \"outputs_identical_across_threads\": %s,\n", all_identical ? "true" : "false");
    std::printf("  \"team_sizes_as_requested\": %s,\n", team_ok ? "true" : "false");
    std::printf("  \"env\": {\"omp_num_procs\": %d, \"openmp_macro\": %d, "
                "\"gmp_version\": \"%s\", \"compiler\": \"%s\", "
                "\"OMP_PROC_BIND\": \"%s\", \"OMP_PLACES\": \"%s\"},\n",
                omp_get_num_procs(), (int)_OPENMP, gmp_version, __VERSION__,
                envs("OMP_PROC_BIND"), envs("OMP_PLACES"));
    std::printf("  \"threads\": [");
    for (size_t i = 0; i < threads.size(); ++i)
        std::printf("%s%d", i ? "," : "", threads[i]);
    std::printf("],\n  \"rows\": [\n");
    bool first = true;
    print_rows(rows_pk, first);
    print_rows(rows_jac, first);
    std::printf("\n  ]\n}\n");

    if (!all_identical) return 3;
    if (!oracle_ok)     return 4;
    if (!team_ok)       return 5;
    return 0;
}
