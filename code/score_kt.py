"""score_kt.py <scored.jsonl> <ends_scored.json> [--out result.json]: the pre-registered kill-test score (PREREG.md).
S = median width over s in [0, 20] mm / median width over s in [25, 40] mm, per end (>= 8 valid samples in each window);
AUC tagged vs control with a 10,000-resample bootstrap over ends (seed 0). Also the same for the peak intensity, and
the AUC per annotator."""
import argparse, collections, json
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument('rows'); ap.add_argument('ends'); ap.add_argument('--out')
a = ap.parse_args()
E = {e['id']: e for e in json.load(open(a.ends))}
by = collections.defaultdict(list)
for ln in open(a.rows):
    r = json.loads(ln)
    if 'width_vox' in r and r['id'] in E: by[r['id']].append(r)


def stat(rows, key):
    s = np.array([r['s_mm'] for r in rows]); v = np.array([r[key] for r in rows])
    n, f = (s >= 0) & (s <= 20), (s >= 25) & (s <= 40)
    if n.sum() < 8 or f.sum() < 8: return None
    return float(np.median(v[n]) / np.median(v[f]))


def auc(pos, neg):
    pos, neg = np.asarray(pos), np.asarray(neg)
    return float(((pos[:, None] > neg[None]).mean() + 0.5 * (pos[:, None] == neg[None]).mean()))


res = dict(per_end={}, dropped=[])
for key in ('width_vox', 'peak'):
    S = {}
    for i, e in E.items():
        st = stat(by.get(i, []), key)
        if st is None:
            if key == 'width_vox': res['dropped'].append(i)
            continue
        S[i] = st
        if key == 'width_vox': res['per_end'][i] = dict(kind=e['kind'], user=e.get('user'), S=round(st, 4))
    K = [S[i] for i in S if E[i]['kind'] == 'kollesis']; C = [S[i] for i in S if E[i]['kind'] == 'control']
    rng = np.random.default_rng(0); bs = []
    for _ in range(10000):
        bs.append(auc(rng.choice(K, len(K)), rng.choice(C, len(C))))
    res[key] = dict(auc=round(auc(K, C), 4), ci95=[round(float(np.percentile(bs, 2.5)), 4), round(float(np.percentile(bs, 97.5)), 4)],
                    n_kollesis=len(K), n_control=len(C), median_S_kollesis=round(float(np.median(K)), 3), median_S_control=round(float(np.median(C)), 3))
    if key == 'width_vox':
        res['per_annotator'] = {u: round(auc([S[i] for i in S if E[i]['kind'] == 'kollesis' and E[i].get('user') == u], C), 4)
                                for u in sorted({E[i].get('user') for i in S if E[i]['kind'] == 'kollesis'})}
a_ = res['width_vox']
res['decision'] = 'GO' if a_['auc'] >= 0.80 and a_['ci95'][0] >= 0.65 else ('WEAK' if a_['auc'] >= 0.65 else 'STOP')
print(json.dumps({k: v for k, v in res.items() if k != 'per_end'}, indent=1))
if a.out: json.dump(res, open(a.out, 'w'), indent=1)
