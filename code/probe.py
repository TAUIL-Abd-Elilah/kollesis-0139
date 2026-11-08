"""probe.py: PHerc 0139 kollesis kill test, CPU only, read-only S3 range reads.

For a fiber end (L0 xyz), walk the fiber's line_points inward from the end (s >= 0, mm) and extend past the end along the
end tangent (s < 0). At each sample, read a 40^3 L1 (4.8 um) crop, take the local sheet normal from the gradient
structure tensor, and measure the papyrus thickness across the sheet: the full width at half height of the intensity
peak nearest the sample, along the normal. Writes one JSON row per sample.

usage: probe.py <ends.json> <out.jsonl> [--step-mm 1] [--in-mm 40] [--out-mm 15] [--threads 24] [--save-xsec DIR]
ends.json: [{"id":..., "file":..., "end": "first"|"last", "kind": "kollesis"|"control"}]
"""
import argparse, json, os, sys, threading, urllib.request
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy import ndimage

B = ('https://vesuvius-challenge-open-data.s3.us-east-1.amazonaws.com/PHerc0139/volumes/'
     '20260102150214-2.399um-0.2m-78keV-masked.zarr/1')
SHAPE = (38477, 13256, 13256); C = 128; UM_L1 = 2 * 2.399; FIB = 'C:/VesuviusTemp/strips/data/human/0139'
_lock = threading.Lock(); _cache = {}


def _get(url, a, b):
    req = urllib.request.Request(url, headers={'Range': f'bytes={a}-{b - 1}'})
    for k in range(4):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception as e:
            if getattr(e, 'code', None) == 404: return None
            if k == 3: raise


def crop(z0, y0, x0, n=40):
    """L1 uint8 crop [z0:z0+n, y0:y0+n, x0:x0+n] via one range request per (chunk, z) plane band."""
    out = np.zeros((n, n, n), np.uint8)
    for cz in range(z0 // C, (z0 + n - 1) // C + 1):
        for cy in range(y0 // C, (y0 + n - 1) // C + 1):
            for cx in range(x0 // C, (x0 + n - 1) // C + 1):
                url = f'{B}/{cz}/{cy}/{cx}'
                za, zb = max(z0, cz * C), min(z0 + n, (cz + 1) * C)
                ya, yb = max(y0, cy * C), min(y0 + n, (cy + 1) * C)
                xa, xb = max(x0, cx * C), min(x0 + n, (cx + 1) * C)
                for z in range(za, zb):
                    off = ((z - cz * C) * C + (ya - cy * C)) * C
                    raw = _get(url, off, off + (yb - ya) * C)
                    if raw is None: continue
                    blk = np.frombuffer(raw, np.uint8).reshape(yb - ya, C)
                    out[z - z0, ya - y0:yb - y0, xa - x0:xb - x0] = blk[:, xa - cx * C:xb - cx * C]
    return out


def measure(p_l1, n=40, half=20):
    """p_l1 = (x, y, z) in L1 voxels. Returns dict with thickness etc."""
    x, y, z = p_l1
    z0, y0, x0 = int(round(z)) - n // 2, int(round(y)) - n // 2, int(round(x)) - n // 2
    v = crop(z0, y0, x0, n).astype(np.float32)
    if (v > 0).mean() < 0.5: return None
    g = ndimage.gaussian_filter(v, 1.0)
    gz, gy, gx = np.gradient(g)
    c = n // 2; sl = (slice(c - 6, c + 7),) * 3
    w = np.ones((13, 13, 13), np.float32)
    J = np.array([[ (a * b)[sl].sum() for b in (gz, gy, gx)] for a in (gz, gy, gx)])
    ev, evec = np.linalg.eigh(J); nrm = evec[:, -1]; anis = (ev[-1] - ev[-2]) / (ev[-1] + 1e-6)
    t = np.arange(-half, half + 0.01, 0.5)
    ctr = np.array([z - z0, y - y0, x - x0])
    pts = ctr[:, None] + nrm[:, None] * t[None]
    prof = ndimage.map_coordinates(g, pts, order=1, mode='nearest')
    # peak nearest t=0 within +-6 vox
    near = np.abs(t) <= 6
    i0 = np.where(near)[0][np.argmax(prof[near])]
    pk = prof[i0]; bg = np.percentile(prof, 10); hh = bg + 0.5 * (pk - bg)
    l = i0
    while l > 0 and prof[l] > hh: l -= 1
    r = i0
    while r < len(t) - 1 and prof[r] > hh: r += 1
    return dict(width_vox=float((r - l) * 0.5), peak=float(pk), bg=float(bg), t_peak=float(t[i0]), anis=float(anis),
                normal=[float(q) for q in nrm], profile=[round(float(q), 1) for q in prof[::2]])


def fiber_path(fname, end):
    d = json.load(open(os.path.join(FIB, fname)))
    P = np.asarray(d['line_points'], np.float64)  # L0 xyz
    if end == 'last': P = P[::-1]
    seg = np.linalg.norm(np.diff(P, axis=0), axis=1); s = np.concatenate([[0], np.cumsum(seg)]) * 2.399e-3  # mm
    return P, s


def samples(fname, end, step, s_in, s_out):
    P, s = fiber_path(fname, end)
    out = []
    for si in np.arange(0, min(s_in, s[-1]) + 1e-9, step):
        out.append((float(si), np.array([np.interp(si, s, P[:, k]) for k in range(3)])))
    k = np.searchsorted(s, 2.0); tang = P[0] - P[min(k, len(P) - 1)]; tang /= np.linalg.norm(tang) + 1e-9
    for so in np.arange(step, s_out + 1e-9, step):
        out.append((-float(so), P[0] + tang * so / 2.399e-3))
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('ends'); ap.add_argument('out')
    ap.add_argument('--step-mm', type=float, default=1.0); ap.add_argument('--in-mm', type=float, default=40)
    ap.add_argument('--out-mm', type=float, default=15); ap.add_argument('--threads', type=int, default=24)
    a = ap.parse_args(); E = json.load(open(a.ends))
    done = set()
    if os.path.exists(a.out):
        for ln in open(a.out): r = json.loads(ln); done.add((r['id'], r['s_mm']))
    jobs = []
    for e in E:
        for s_mm, p in samples(e['file'], e['end'], a.step_mm, a.in_mm, a.out_mm):
            if (e['id'], round(s_mm, 3)) not in done: jobs.append((e, round(s_mm, 3), p))
    print(len(E), 'ends', len(jobs), 'samples to do', flush=True)
    fo = open(a.out, 'a')

    def run(j):
        e, s_mm, p = j
        try:
            m = measure(p / 2.0)
        except Exception as ex:
            m = dict(error=str(ex)[:120])
        row = dict(id=e['id'], kind=e['kind'], s_mm=s_mm, l0_xyz=[round(float(q), 1) for q in p], **(m or dict(empty=True)))
        with _lock: fo.write(json.dumps(row) + '\n'); fo.flush()
    with ThreadPoolExecutor(a.threads) as ex:
        for i, _ in enumerate(ex.map(run, jobs)):
            if i % 200 == 0: print(i, flush=True)
    print('done', flush=True)
