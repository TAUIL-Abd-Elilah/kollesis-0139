# kollesis-0139

Can kolleseis, the joins where papyrus sheets were pasted together, be found in the CT and used as full-height same-winding constraints for the spiral fit? The Vesuvius Challenge winding-annotation page lists this under "Things we have not tried".

**Status, 8 Nov 2026: kill test pre-registered, run in progress.** `PREREG.md` fixes the data, measurement, statistic and decision rule before the scored run. The first step asks whether the CT shows a papyrus-thickness signature at the 78 joins tagged by the PHerc. 0139 annotators, compared with untagged fiber ends. A detector and a fit A/B follow only if it passes.

- `code/probe.py`: thickness along fibers near their ends, from public S3 L1 data, using range reads only.
- `data/ends_scored.json`, `data/ends_pilot.json`: the fiber ends (tagged and control) and the pilot ends excluded from scoring.

Data: PHerc. 0139 volume 20260102150214, from the Vesuvius Challenge open data. Fibers and `kollesis_termination` tags come from the 0139 annotation team (lt, kb, dj): https://dl.ash2txt.org/other/annotators/0139_fibers/

Credits: the PHerc. 0139 annotation team; ChaseHendrick (#2386, reading kollesis tags); skyward7187 (#2002/#2185, kollesis seams in VC3D); the Vesuvius Challenge team (the kolleisis idea in the winding-annotation docs).

Code MIT. Developed with AI assistance under my direction.
