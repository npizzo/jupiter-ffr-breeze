# Direct day-1-to-50 continuations

At the user's explicit request, remove the intermediate day-2 timing gate and run the three saved 256 x 256 x 64 cases directly to day 50 on A100/H100 GPUs. The staged continuation jobs had not started and consumed zero GPU time; they were held before replacement. No completed simulation interval is repeated.

All physical and numerical configuration values remain unchanged. Each job verifies original configuration/source integrity, requires A100/H100, loads its retained day-1 checkpoint and sets stop_time to 50 Earth days. There is no forecast-based stop. Ordinary CFL, NaN and cap checks remain, as do periodic checkpoints, three-hour movie output, daily full fields, and the 600-second heat controller.

The allocation remains 5h50m per continuation. Adding the prior 568/564/553 seconds stays below the original six GPU-hours per case and eighteen total. The normal internal wall limit leaves five minutes to save state; it is not an intermediate timing test. No simulation source or prior output is modified.

The new root links the original day0_1 directory and writes a single day1_50 directory. The dependent renderer reuses the existing two-segment timestamp, restart-frame, shared-color-scale and full-decode checks. Compute and rendering remain on Unity. The original staged jobs and their pending renderer are replaced, not run in addition.
