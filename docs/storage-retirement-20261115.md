# Historical output retirement, 15 November 2026

The user authorized removing older runs to recover shared Unity work storage. An audited manifest retired 680 files totaling 87,128,849,101 bytes. Reported free capacity subsequently increased from about 111 GiB to 192 GiB in the nominal 1 TB work allocation. The allocation also contains work outside this repository; this campaign occupied about 365 GiB before cleanup.

Only middle `frame_*.jld2` files from the completed, published `large40-v1` and `shear20-v2` campaigns, and JLD2 files from the superseded `fplane-point4-v1/preflight`, were removed. All historical production checkpoints, first/last frames of each affected segment, CSV/TOML/JSON diagnostics, published movies and figures were retained. Current f-plane and beta-plane campaigns, their checkpoints and the active refinement inputs were untouched.

The older gallery assets were previously verified under publication commit `e976812f77878357469b89053aa48b7ab1b73a5a`. Those published movies remain available; the retired intermediate raw frames can no longer be reanalyzed directly. Reproduction requires the retained configurations/checkpoints and rerunning the relevant interval.

The exact retirement plan and receipt remain on Unity under `docs/cleanup-20261115-plan.json` and `docs/cleanup-20261115-receipt.json`, including paths, sizes and pre-deletion modification times. No output was downloaded for cleanup.
