# Test suite slowness and memory use (FTA)

Target commit: 3b51672 (root fillet included). Investigated by an engineer agent in a git worktree; this file is the main agent's transcription of its report.


Scope: the full test suite of wave-ball-reducer (`pixi run test`, 19 tests) is slow (761 s before the root fillet was added to the gear teeth, 1257 s after) and peaks at about 20.7 GB of resident memory. All measurements below were taken at the worktree HEAD (commit 3b51672, which already includes the fillet), with the main checkout's pixi Python. Every Python process ran inside `systemd-run --user --scope -p MemoryMax=24G -p MemorySwapMax=0`, one heavy process at a time. None hit the cap.

#### Measurement 1: per-test durations (`pytest --durations=0`, run in full, 19 passed)

| Test | Time [s] |
| --- | --- |
| test_planetary: no interference and tooth gap, pose 1 (planet 1, input 0.9 rad) | 336.3 |
| same, pose 2 (planet 2, input 2.3 rad) | 327.8 |
| same, pose 0 (planet 0, input 0) | 326.3 |
| test_reducer: full input turn (12 poses) | 144.3 |
| test_planetary: build and reload | 81.6 |
| test_reducer: rebuild with nb=15 | 25.7 |
| test_planetary: module fixture `home` (make_parts, setup) | 22.0 |
| test_reducer: build and reload | 8.0 |
| test_reducer: fixture `home` (setup) | 3.3 |
| all others (11 tests) | under 1.3 each |
| **Total** | **1279.4** (peak RSS 20.6 GB) |

The three planetary mesh tests take 990 s, 77 % of the total.

#### Measurement 2: breakdown inside one planetary pose (probe script, input angle 0.9 rad)

| Step | With fillet (current) | Fillet disabled (`_fillet` made identity) | Ratio |
| --- | --- | --- | --- |
| Faces of one planet / sun / housing / ring_lid | 1253 / 633 / 1492 / 1461 | 869 / 441 / 1012 / 981 | about 1.44x |
| make_parts (all 12 parts) | 21.8 s, RSS 8.7 GB | 13.1 s, RSS 5.8 GB | 1.7x |
| interference, all 66 part pairs | 52.3 s | 29.4 s | 1.8x |
| of which housing vs ring_lid (both fixed, same result in every pose) | 20.6 s | 10.0 s | 2.1x |
| of which each planet vs housing / vs ring_lid | about 2.9 s each | about 1.7 s each | 1.7x |
| BREP distance sun-planet | 80.2 s | 45.6 s | 1.8x |
| BREP distance housing-planet | 101.3 s | 53.9 s | 1.9x |
| BREP distance ring_lid-planet | 98.0 s | 52.7 s | 1.9x |
| Peak RSS after interference | 19.7 GB | 12.0 GB | 1.6x |
| `build` (make_parts plus STEP/STL export) | 28.5 s | not measured | |

One pose test = interference 52 s + three distances 280 s = about 330 s, which matches the 326 to 336 s in measurement 1. Across the three poses the nine BREP distance calls cost about 840 s (66 % of the whole suite) and interference about 157 s.

#### Measurement 3: steel-ball reducer (same probe)

make_parts 3.0 s (RSS 0.66 GB); faces per part 4 to 67; interference over 55 pairs in one pose 7.7 s; build 5.1 s; peak RSS 1.1 GB. The 144 s of the full-turn test is 12 poses of about 8 s interference (about 92 s) plus ball-to-disc distance checks and the 400,000-point groove arrays. Its parts are simple, so it is unrelated to the fillet regression.

#### Measurement 4: memory

| Probe | Result |
| --- | --- |
| One planet gear solid (`herringbone`, 624-point outline) built alone | 5.2 s, RSS grows by 8.37 GB |
| Ring cutter built right after it | RSS grows by 0.00 GB (reuses the freed heap) |
| After deleting both solids | RSS stays at 8.72 GB |
| After `malloc_trim(0)` from libc | RSS drops to 0.45 GB |
| `del` of all parts and `gc.collect()` after the interference loop | RSS stays at 19.68 GB |

So the memory is mostly transient working memory of the twisted extrusion and of the booleans. glibc keeps the freed heap and does not return it to the OS. The process peak (about 20 GB) is reached during the planetary interference loop; the largest single boolean is housing vs ring_lid. Live data (the 12 solids themselves) is small. Holding all Solids is not the cause. I did not separate how much of the 8.4 GB comes from `Solid.extrude_linear_with_rotation` and how much from the following `fuse` of the two halves; that remains unconfirmed.

#### Root cause

1. Main cost: the gear flanks are a polyline of 624 points per planet. Each segment becomes its own twisted (ruled B-spline) face, so a planet has about 1250 faces and housing and ring_lid about 1500. OCCT's BREP distance (`distance_to`) and boolean common scale super-linearly with face count. The test calls distance 9 times and booleans 66 times per pose on these solids.
2. Regression mechanism: the fillet replaces each root corner with a 5-point arc. This raises the face count by 44 % and roughly doubles the distance and boolean times (1.8 to 1.9x measured). The 761 s to 1257 s jump matches: nine distances +about 390 s, interference +about 70 s, part builds +about 30 s.
3. Memory: the twisted extrusion and booleans of these many-face solids take GB-scale temporary allocations (8.4 GB for one planet). glibc keeps the freed memory, so RSS climbs to about 20 GB and stays there. With pytest-xdist each worker would carry this, so 3 workers would need about 60 GB.

#### Proposals

| Option | Estimated saving | Coverage loss | Effort |
| --- | --- | --- | --- |
| A. Measure flank clearance on 2D sections (for example the `section()` at 3 to 5 heights through the gear) instead of 3D BREP distance on the full helical solids. Per the `herringbone` docstring, slices at the same height are a valid mesh. | about 800 s (840 s now; 2D distances take seconds) | Checked only at the sampled heights. Clearance between sampled heights is not checked. | small (a helper and the test loop) |
| B. `slow` marker, excluded by default (`-m "not slow"` in pytest.ini) on the 3 planetary mesh tests and the reducer full-turn test; run them explicitly before commits | about 1135 s off the default run (remaining run about 145 s) | None if the slow set is still run before commits; otherwise the mesh checks run less often | very small |
| C. Remove waste: (1) check housing vs ring_lid once, since both are fixed and the result is the same in every pose (saves about 40 s of the 60 s); (2) in test_planetary build-and-reload, compare against `parts.json` or the `home` fixture instead of calling `make_parts(Params())` a second time (saves about 22 s) | about 60 to 65 s | none | very small |
| D. Coarser tooth profile for tests (smaller `pts_per_flank` from 8, fewer arc points in `_fillet` from 4 segments) or a fillet expressed as fewer points | proportional to face count; going back to the pre-fillet count gives about 45 % of the gear-related time (about 500 s) | Tests no longer check the exact shipped geometry | small |
| E. Memory: call `malloc_trim(0)` after building parts or between tests (for example in an autouse fixture), or set `MALLOC_TRIM_THRESHOLD_` | RSS between tests about 20 GB down to under 1 GB (measured 8.7 to 0.45 GB); peak during the heaviest boolean unchanged | none | very small |
| F. pytest-xdist (estimate only, not measured as instructed): 3 workers would run the 3 pose tests in parallel, about 1279 to about 560 s | about 700 s | none | small, but each worker peaks near 20 GB, so OOM risk on 62 GB. Unsafe without A and E. |
| G. Bounding-box pre-filter in `interference` (XY as well as Z) | small: gear pairs overlap in their bounding boxes, so almost none are skipped | none | small |
| H. Approximate distance (mesh or point sampling) | close to A | tolerance depends on the tessellation; for a 0.3 mm clearance limit it needs a fine mesh | medium |

Top 3 recommendations:
1. A (2D section distance): about 800 s saved, the largest single gain, small coverage loss.
2. C (pose-invariant pair once, no duplicate make_parts): about 60 s, no coverage loss.
3. B (slow marker) for everyday runs: default run about 145 s. Or E if memory is the concern; it makes later parallel runs safer.
A and C together: about 1279 s to about 420 s for the full suite with nearly the same coverage. Adding B makes the default run about 145 s.

#### Fault tree analysis log (full text)

Phase 0 triage:
- Linter or type checker catches it: no (a performance issue, not a static error).
- git log shows a culprit: partially. Commit d102fed added the planetary reducer, and the fillet is in HEAD. A cause is suggested but not a mechanism, so the answer is no.
- Self-explanatory error: no (no error, only slowness).
Result: full tree.

Phase 1 top event: `pytest tests` takes 1279 s (reproduced here; the reported value was 1257 s) and peaks at 20.6 GB RSS, 19 of 19 passing. Reproducible: a single full run.

Phase 2 tree (OR gate; hypotheses written from structure only):
- H1 Logic or duplicate work: parts rebuilt repeatedly (fixtures, build, rebuild).
- H2 Data and geometry: the many-face flank solids make booleans (interference) slow; the fillet raised the face count.
- H3 Logic and algorithm: the BREP distance calls are slow on helical solids.
- H4 I/O: STEP/STL export and re-import are slow.
- H5 Environment: CPU contention with another test run.
- H6 Dependency: an OCCT or build123d version change.
- H7 Race condition: excluded. The tests are single-process with no concurrency.
- H8 Resource and state: memory retention (all solids kept, boolean results piling up, or allocator retention) makes it slow or memory-heavy.
- Configuration: excluded. `pytest.ini` only sets pythonpath, and no markers or flags are involved.

Phase 3 order:
1. Per-test durations; (2) cheap and (3) high information value, since it splits H1/H2/H3/H4 at once.
2. Per-step probe with fillet and without fillet; (3), since it isolates H2/H3 and tests the fillet regression.
3. RSS at each step; (3) for H8.
4. H5 and H6 from the environment; (2) cheap.

Phase 4:
- H1: predicted that if true, setup or rebuild time would be a large share. Observed: planetary `home` setup 22 s once (module scope, shared by the 3 pose tests), build-and-reload 82 s (of which make_parts runs twice, about 44 s), reducer rebuild 26 s. About 100 s in all, 8 % of 1279 s. Contributing, not the root cause.
- H2: predicted that interference time would scale with the face count. Observed: 52.3 s with fillet vs 29.4 s without, for faces 1.44x. Interference is about 157 s in total. Contributing.
- H3: predicted that the distance calls would dominate if true. Observed: 80 to 101 s per call with fillet, 46 to 54 s without; 9 calls give about 840 s = 66 % of the total. Confirmed as the main cost.
- H4: predicted that export would be a large share. Observed: planetary build 28.5 s including make_parts 22 s, so export about 7 s; reducer build 5 s. Eliminated as a main cause.
- H5: no other pytest process was running during the final measurement (checked with `pgrep`; 46 GB free). The probe's sum per pose (332 s) equals the pytest time per pose (326 to 336 s). Eliminated.
- H6: the same environment was used for the with-fillet and without-fillet runs, and the ratio follows the face count. The dependency cannot explain the fillet delta. Eliminated as the regression cause, though not checked against older versions.
- H8: predicted that if holding the solids were the cause, RSS would drop after `del`. Observed: no drop after `del`, but a drop to 0.45 GB after `malloc_trim`; one planet alone allocated 8.37 GB transiently. Confirmed as allocator retention of transient peaks. Holding solids is eliminated.

IS / IS-NOT:
- IS: planetary pose tests with fillet (about 330 s, 20 GB).
- IS-NOT 1: the same steps with the fillet disabled (about 196 s per pose equivalent, 12 GB).
- IS-NOT 2: the reducer (faces at most 67): 8 s per pose, 1.1 GB.
- Face count is present in every IS case and smaller in every IS-NOT case, so it explains both columns.
- I/O (H4) is present in both, so it is eliminated.

Phase 5: Reproduction (full run), Isolation (`herringbone` polyline faces + `distance_to` / boolean in the planetary pose test), Mechanism (faces 1.44x gives time 1.8x; transient GB allocations kept by glibc), Exclusion (above). Unresolved: the split of the 8.4 GB between the twisted extrusion and the `fuse` in `herringbone`.

Cleanup: no repository file was changed (`git status` clean); `grep -rn '\[FTA H' .` in the worktree returns nothing. The probe scripts are only in the scratchpad.


## Follow-up (2026-11-12): better countermeasures

Probe: sun + planet only (626 / 1250 faces), pose input angle 0.9 rad, sun–planet pair, `BRepExtrema_DistShapeShape` directly. Run serially under `systemd-run --user --scope -p MemoryMax=24G -p MemorySwapMax=0`, 20-core machine.

| mode | distance [mm] | time [s] | peak RSS [GB] |
|---|---|---|---|
| base (= `distance_to`) | 0.094717 | 80.2 | 9.0 |
| `SetMultiThread(True)` | 0.094717 | 10.0 | 9.0 |
| face filter (faces within partner tip radius + 0.5 mm; 134/626, 146/1250) | 0.094717 | 9.6 (+0.13 filter) | 8.9 |
| filter + multithread | 0.094717 | 2.1 | 9.0 |

Peak RSS comes from building the planet solid, not from the distance call. build123d 0.13.0 `distance_to` does not expose `SetMultiThread`; booleans already use `SetRunParallel(True)`.

Conclusion: multithread + tip-circle face filter replaces option A. It is exact while the true gap < 0.5 mm margin (the test already requires < 0.3 mm). Estimated 9 distance calls go from ~840 s to ~20–30 s. Estimated suite with 1+2 + C + E: ~400 s, vs ~420 s for A + C + E, and no loss of coverage between slices. Not yet measured: housing–planet and ring_lid–planet pairs (they need a different filter shape). Also considered: skipping booleans for pairs already proven gap > 0 (~−27 s); spline flanks (rejected, product geometry change, net effect unknown); 2D sweep of the outline (checks the profile function, not the solid); mesh distance (new dependency); part caching (already effectively done by module-scoped `home`).

Source: https://dev.opencascade.org/doc/refman/html/class_b_rep_extrema___dist_shape_shape.html

## Result (2026-11-12): fixes applied

Applied: multithread + face filter (`gear_gap` in tests/test_planetary.py), fixed-part pairs checked once (`test_fixed_parts`, `interference(skip=...)` in reducer.py), `home` fixture reused in `test_build_and_reload`, `malloc_trim(0)` after `make_parts`. planetary.py untouched.

- Ring pair probe (housing–planet_1, pose 0.9): plain `distance_to` 0.09014180302024451 mm / 96.57 s; multithread 20.98 s; filter + multithread 4.01 s, identical distance.
- Full suite under 24G cap: 20 passed, 403.07 s (was 19 passed, 1264.64 s). Peak RSS 19.7 GB (was 20.7 GB): the peak comes from building/booleans, not distance.
- test_no_interference_and_tooth_gap: ~40 s per pose (+21 s fixture). Slowest remaining: test_reducer::test_full_input_turn 141.6 s, test_build_and_reload 59.0 s.
