# Resume the three saved 256-squared cases on A100/H100

The user explicitly authorized this continuation after the original L40S jobs completed one day but failed the runtime-budget gate. The saved day-1 checkpoints and output are retained; day 0-1 is not rerun.

- Same 256 x 256 x 64 grid, 6000-km box, 1200-km eddies, random RMS 5.5022 m/s, random Ro=0.1, no/cyclonic/anticyclonic initial shear, heating, cap, boundaries, viscosities, numerical methods and output/controller cadence as `strong256_v1`.
- Only hardware eligibility and restart staging change: A100/H100 only, advance day 1 to day 2 to measure performance, then continue the day-2 checkpoint to day 50 if the cost gate passes. No new initial conditions.
- Three independent GPU jobs. Prior allocations were 568, 564, and 553 GPU-seconds. Each continuation requests at most 21000 seconds (5h50m); prior plus new allocation is less than six GPU-hours per case, preserving the eighteen-GPU-hour combined ceiling. The internal continuation includes a 25% runtime margin and stops five minutes before the allocation limit.
- Validate original configuration and frozen source hashes, CUDA runtime, GPU model, checkpoint hash, restart time, and saved restart fields. Existing driver checks finite state, CFL and cap stability. These checks do not prove spatial convergence or late-time fidelity.
- Failed speed checks retain day 2 and do not retry automatically. Successful runs produce separate day1_2 and day2_50 directories; the old day0_1 remains in its original location.
- A replacement dependent Unity CPU renderer joins the three segments, checks both duplicate restart frames and all 401 three-hourly timestamps, renders the three-depth vorticity/tracer movies with shared fixed scales, and fully decodes them for validation. The original dependent renderer was cancelled after its predecessor gates failed; it is not reused.

No computation or rendering runs on the local machine. Model source and previous campaign files remain unchanged.
