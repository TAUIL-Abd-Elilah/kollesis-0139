# Kollesis kill test, PHerc 0139 (pre-registered 8 Nov 2026, after a 12-end pilot, before the scored run)

**Question.** Does the CT show a papyrus-thickness signature at the joins that the 0139 annotators tagged, strongly enough to justify building a kollesis detector? (Read-only; CPU; S3 range reads of L1, 4.8 um.)

**Data.**
- 78 `kollesis_termination` fiber ends in the public 0139 annotator fibers (lt 35, kb 31, dj 12; 77 H fibers).
- Controls: 78 untagged ends of H fibers, drawn at random (seed 0) from 430. Ends tagged `break` are excluded.
- The 12 pilot ends (6 + 6, `ends_pilot.json`) are excluded from scoring. Scored set: `ends_scored.json`, 72 + 72.

**Measurement** (`probe.py`, frozen). Samples run every 1 mm along the fiber, from the end inward to 40 mm, and 15 mm past the end along the end tangent. At each sample, the papyrus thickness is the half-height width of the intensity peak nearest the sample, along the structure-tensor normal (L1 voxels).

**Statistic per end.** S = median width over s in [0, 20] mm divided by median width over s in [25, 40] mm. An end is scored only if it has at least 8 valid samples in each window; ends dropped this way are reported.

The pilot motivated this statistic. The 6 pilot tagged ends had S of about 1.2-1.7; the 6 pilot controls had S of 0.5-1.5.

**Primary measure.** AUC of S, tagged ends vs control ends, with a 95 % bootstrap CI over ends (10,000 resamples, seed 0).

**Decision.**
- **GO** (build a detector and run a constraint fit A/B): AUC >= 0.80 and CI lower bound >= 0.65.
- **WEAK** (report it; try one other feature only if it is pre-registered first): 0.65 <= AUC < 0.80.
- **STOP** (record the null): AUC < 0.65.

**Also reported:** the same AUC for the intensity peak instead of width, per annotator, and per-end S values.
482f4c9885121cb41419413b87e79444047a0de337a80ff4a533b73b921bfe95 *probe.py
435d3e031ae7c82d07730acd8ba77fe703439791ce0f7d8db84c684da5573acc *ends_scored.json
9d4c309721d6d98278f417c3535b21331391d62811bd8408b5ca589e5b22f9e1 *ends_pilot.json
