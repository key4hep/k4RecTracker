#!/usr/bin/env python3
# Compare two reconstructed track samples against generator-level truth.
#
# Usage: compare_tracking.py [options] <fileA.root[:collectionA]> <fileB.root[:collectionB]>
#
# Defaults: A = "FittedTracks" (GGTF + Genfit), B = "SiTracksCT_Refitted"
# (conformal tracking + RefitFinal).
# Truth selection (--truth):
#   gun    - generator mu- from the muon gun sample
#   stable - prompt stable charged particles pT > 1 GeV (for physics events, e.g. Z->uds)
#
# Matching: a track is matched to a truth particle if |dphi| < 0.02 rad and
# |d(tanLambda)| < 0.02 using the track state at the IP (greedy best match,
# each track matches at most one truth particle). pT is derived from omega
# assuming the nominal IDEA 2 T solenoid field.

import argparse
import json
import math
import statistics

from podio.reading import get_reader

MATCH_DPHI = 0.02
MATCH_DLAMBDA = 0.02
BFIELD_TESLA = 2.0  # nominal IDEA solenoid field
PT_PER_OMEGA = 0.299792458e-3 * BFIELD_TESLA  # GeV for pT = c*B/|omega|, omega in 1/mm

AT_IP = 1  # edm4hep::TrackState::AtIP

DIGI_COLLECTIONS = ["VTXBDigis", "VTXDDigis", "SiWrBDigis", "SiWrDDigis"]

TRUTH_MODE = "gun"  # or "stable"


def select_truth(frame):
    """Generator-level truth particles for the efficiency denominator."""
    truth = []
    for p in frame.get("MCParticles"):
        if p.getGeneratorStatus() != 1 or p.getCharge() == 0:
            continue
        m = p.getMomentum()
        pt = math.hypot(m.x, m.y)
        if pt < 1.0:
            continue
        v = p.getVertex()
        rho = math.hypot(v.x, v.y)
        if TRUTH_MODE == "gun":
            if abs(p.getPDG()) != 13 or rho > 1.0:  # only the prompt gun muons
                continue
        else:  # stable: any prompt stable charged particle (beam smear is um-scale)
            if rho > 0.5:
                continue
        truth.append(
            dict(
                pt=pt,
                p=math.sqrt(m.x**2 + m.y**2 + m.z**2),
                phi=math.atan2(m.y, m.x),
                lam=m.z / pt,  # tanLambda at the origin
                matched=False,
            )
        )
    return truth


def state_at_ip(track):
    states = track.getTrackStates()
    for s in states:
        if int(s.location) == AT_IP:
            return s
    return states[0] if len(states) else None


def wrap_phi(dphi):
    while dphi > math.pi:
        dphi -= 2 * math.pi
    while dphi < -math.pi:
        dphi += 2 * math.pi
    return dphi


def analyze(path, track_collection):
    reader = get_reader(path)
    events = reader.get("events")
    assert len(events) > 0, f"{path}: no events"

    out = {"file": path, "collection": track_collection, "events": [], "residuals": {}}
    counters = dict(
        events=0,
        muons=0,
        matched=0,
        tracks=0,
        tracks_matched=0,
        failed=0,
        no_ip_state=0,
        duplicates=0,
    )
    res = {k: [] for k in ["rel_pt", "d0", "z0", "phi", "tanLambda", "hits", "chi2ndf"]}
    hits_by_collection = {c: 0 for c in DIGI_COLLECTIONS}

    for iev, frame in enumerate(events):
        counters["events"] += 1
        digi_ids = {}
        for c in DIGI_COLLECTIONS:
            try:
                for h in frame.get(c):
                    digi_ids[(h.id().collectionID, h.id().index)] = c
            except Exception:
                pass

        muons = select_truth(frame)
        counters["muons"] += len(muons)

        tracks_ev = []
        for t in frame.get(track_collection):
            counters["tracks"] += 1
            nhits = len(t.getTrackerHits())
            res["hits"].append(nhits)
            for h in t.getTrackerHits():
                c = digi_ids.get((h.id().collectionID, h.id().index))
                if c:
                    hits_by_collection[c] += 1
            if t.getChi2() < 0 or t.getNdf() <= 0:
                counters["failed"] += 1
                tracks_ev.append(None)
                continue
            ndf = t.getNdf()
            res["chi2ndf"].append(t.getChi2() / max(ndf, 1))
            s = state_at_ip(t)
            if s is None or (s.location != AT_IP and abs(s.referencePoint[0]) > 100):
                counters["no_ip_state"] += 1
            info = None
            if s is not None:
                omega = s.omega
                info = dict(
                    phi=s.phi,
                    lam=s.tanLambda,
                    pt=PT_PER_OMEGA / abs(omega) if omega else 0.0,
                    omega=omega,
                    d0=s.D0,
                    z0=s.Z0,
                    at_ip=(int(s.location) == AT_IP),
                )
            tracks_ev.append(info)

        # greedy matching by normalized distance
        candidates = []
        for im, mu in enumerate(muons):
            phi_mc = mu["phi"]
            for it, info in enumerate(tracks_ev):
                if info is None:
                    continue
                dphi = wrap_phi(info["phi"] - phi_mc)
                dlam = info["lam"] - mu["lam"]
                if abs(dphi) < MATCH_DPHI and abs(dlam) < MATCH_DLAMBDA:
                    candidates.append((dphi**2 + dlam**2, im, it, dphi, dlam))
        candidates.sort()
        used_mc, used_trk = set(), set()
        ev_match = 0
        for d, im, it, dphi, dlam in candidates:
            if im in used_mc or it in used_trk:
                if im in used_mc and it not in used_trk:
                    counters["duplicates"] += 1
                continue
            used_mc.add(im)
            used_trk.add(it)
            ev_match += 1
            mu = muons[im]
            info = tracks_ev[it]
            res["rel_pt"].append((info["pt"] - mu["pt"]) / mu["pt"])
            res["phi"].append(dphi)
            res["tanLambda"].append(dlam)
            if info["at_ip"]:
                res["d0"].append(info["d0"])
                res["z0"].append(info["z0"])
        counters["matched"] += ev_match
        counters["tracks_matched"] += len(used_trk)
        out["events"].append(
            dict(
                event=iev,
                muons=len(muons),
                tracks=sum(1 for x in tracks_ev if x is not None),
                failed=sum(1 for x in tracks_ev if x is None),
                matched=ev_match,
            )
        )

    def stats(v):
        if not v:
            return {"n": 0}
        d = {"n": len(v), "mean": statistics.fmean(v)}
        if len(v) > 1:
            d["rms"] = statistics.stdev(v)
        return d

    counters["unmatched_tracks"] = (
        counters["tracks"] - counters["failed"] - counters["tracks_matched"]
    )
    out["counters"] = counters
    for k, v in res.items():
        out["residuals"][k] = stats(v)
    out["raw"] = res  # per-matched-track values, used by plot_tracking.py
    total_hits = sum(hits_by_collection.values()) or 1
    out["hit_fraction"] = {c: n / total_hits for c, n in hits_by_collection.items()}
    out["residuals"]["hits"]["mean_per_track"] = (
        statistics.fmean(res["hits"]) if res["hits"] else 0
    )
    return out


def fmt(tag, r, scale=1.0, unit=""):
    if r.get("n", 0) == 0:
        return f"{tag:>8s}: n=0"
    return f"{tag:>8s}: n={r['n']:3d}  mean={r.get('mean', 0) * scale:+.4f}  rms={r.get('rms', 0) * scale:.4f} {unit}"


def report(tag, out):
    c = out["counters"]
    eff = c["matched"] / c["muons"] if c["muons"] else 0
    eff_err = math.sqrt(eff * (1 - eff) / c["muons"]) if c["muons"] else 0
    r = out["residuals"]
    print(f"\n=== {tag}  ({out['file']} :: {out['collection']})")
    print(f"events={c['events']} truth particles={c['muons']} (mode={TRUTH_MODE})")
    print(
        f"tracks total={c['tracks']} failed fits={c['failed']} matched={c['tracks_matched']} "
        f"unmatched={c['unmatched_tracks']} dup-candidates={c['duplicates']}"
    )
    print(f"truth matching efficiency = {eff:.3f} +/- {eff_err:.3f}")
    print(
        f"hits/track mean = {r['hits']['mean_per_track']:.1f};  hit content: "
        + ", ".join(f"{k}={v:.0%}" for k, v in out["hit_fraction"].items())
    )
    print(fmt("rel_pt", r["rel_pt"]))
    print(fmt("d phi", r["phi"], 1e3, "mrad"))
    print(fmt("d tanL", r["tanLambda"], 1e3))
    print(fmt("d0", r["d0"], 1, "mm"))
    print(fmt("z0", r["z0"], 1, "mm"))
    print(fmt("chi2/ndf", r["chi2ndf"]))
    return eff


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        "A", help="fileA.root or fileA.root:collectionA (default collection FittedTracks)"
    )
    ap.add_argument(
        "B", help="fileB.root or fileB.root:collectionB (default collection SiTracksCT_Refitted)"
    )
    ap.add_argument("--json", default=None, help="output json file (default compare_results.json)")
    ap.add_argument(
        "--truth", choices=["gun", "stable"], default="gun", help="truth particle selection"
    )
    args = ap.parse_args()

    global TRUTH_MODE
    TRUTH_MODE = args.truth

    def split(spec, default_col):
        if ":" in spec:
            return spec.split(":", 1)
        return spec, default_col

    fa, ca = split(args.A, "FittedTracks")
    fb, cb = split(args.B, "SiTracksCT_Refitted")
    outa = analyze(fa, ca)
    outb = analyze(fb, cb)
    report("Pipeline A: GGTF + Genfit", outa)
    report("Pipeline B: ConformalTracking + Refit", outb)
    json_path = args.json or f"compare_results_{TRUTH_MODE}.json"
    with open(json_path, "w") as f:
        json.dump({"A": outa, "B": outb}, f, indent=2)
    print(f"\nPer-event detail and per-track values are in {json_path}")
    for tag, out in (("A", outa), ("B", outb)):
        evs = out["events"]
        flat = " ".join(f"{e['matched']}/{e['muons']}({e['tracks']}+{e['failed']}F)" for e in evs)
        print(f"per-event matched/truth(tracks) [{tag}]: {flat}")


if __name__ == "__main__":
    main()
