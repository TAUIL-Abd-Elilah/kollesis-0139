# kollesis-0139

Can kolleseis, the joins where papyrus sheets were pasted together, be found in the CT and used as full-height same-winding constraints for the spiral fit? The Vesuvius Challenge winding-annotation page lists this under "Things we have not tried".

**Status, 9 Nov 2026: kill test scored. Verdict: STOP (pre-registered null).** The rule in `PREREG.md` was fixed and pushed before the scored run. With it, the thickness signature does not separate the tagged joins from untagged fiber ends:
- AUC 0.538, 95 % CI 0.427-0.651, on 60 tagged ends vs 42 controls. The ends dropped for too few samples are listed in `results/kill_test.json`.
- Median S is 1.04 for tagged ends and 0.97 for controls.
- The peak-intensity variant gives AUC 0.542. Per annotator: dj 0.59, kb 0.60, lt 0.47.

The 12-end pilot looked strong (6 of 6 tagged ends had S >= 1.2). It did not hold on the scored set, which is why the pilot was excluded and the bar was set in advance.

What this does and does not show: a simple half-height thickness, measured along the fiber within 20 mm of a tagged end at 4.8 um, carries no usable join signal here. It does not rule out a learned detector, other features, other resolutions, or joins located some other way. Following the pre-registration, no other feature was tried on these data.

- `code/probe.py`: thickness along fibers near their ends, from public S3 L1 data, using range reads only.
- `code/score_kt.py`: the pre-registered score.
- `results/kill_test.json`: the full result, including S per end.
- `data/ends_scored.json`, `data/ends_pilot.json`: the fiber ends (tagged and control) and the pilot ends excluded from scoring.

Data: PHerc. 0139 volume 20260102150214, from the Vesuvius Challenge open data. Fibers and `kollesis_termination` tags come from the 0139 annotation team (lt, kb, dj): https://dl.ash2txt.org/other/annotators/0139_fibers/

Credits: the PHerc. 0139 annotation team; ChaseHendrick (#2386, reading kollesis tags); skyward7187 (#2002/#2185, kollesis seams in VC3D); the Vesuvius Challenge team (the kolleisis idea in the winding-annotation docs).

Code MIT. Developed with AI assistance under my direction.
