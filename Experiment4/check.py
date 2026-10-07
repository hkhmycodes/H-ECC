import json, glob, sys

files = sys.argv[1:] if len(sys.argv) > 1 else sorted(glob.glob("results/exp4_*.json"))
if not files:
    print("no files matched")
    sys.exit(0)

for f in files:
    try:
        d = json.load(open(f))
    except Exception as e:
        print(f"{f}: {e}")
        continue
    print(f"\n=== {f} ===")
    print(f"n={d['n']}  seed={d['seed']}  repeats={d['repeats']}")
    for r in d["rows"]:
        s = r["samples_s"]
        if len(s) < 2:
            continue
        spread = (max(s) - min(s)) / min(s) * 100
        flag = "  <-- noisy" if spread > 20 else ""
        print(f"  {r['kind']:4s} T={r['T']:2d}  "
              f"median={r['T_median_s']:8.3f}s  "
              f"speedup={r['speedup']:.3f}x  "
              f"eff={r['efficiency']:.3f}  "
              f"spread={spread:5.1f}%{flag}")
