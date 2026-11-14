# Recovery of the f-plane analysis environment

The v2 GPU preflight (1725453) completed all eight simulation pilots and passed the revised runtime and conservative-remapping checks. The job then failed when its Python diagnostics tried to import SciPy from the shared gallery environment, where SciPy was absent. Its production and review dependencies were automatically cancelled before any production began.

This repair leaves all v2 physical source, configurations, checkpoints and completed GPU pilots unchanged. A separate Unity Python environment, `tools/point4-analysis-env`, adds SciPy. Its full package lock is [recorded here](../../environments/fplane-point4-analysis-v3.txt); the original gallery environment is preserved.

Recovery job 1737737 verifies the existing GPU gate, source/configuration hashes and restart-bundle hash, then runs the previously unreachable analytic and saved-field diagnostics. All 13 checks passed. Only then is the combined gate written. A separate end-to-end smoke job (1737812) exercises the full analysis, plots, three-depth rendering and complete movie decode on the actual heated pilot. The resubmission launcher requires that movie check as well before scheduling production.

The new production array still uses the frozen v2 GPU driver and 16-hour allocation, with at most two GPUs running simultaneously. Review and comparison jobs use the separate, complete Python environment. Failed attempts remain available; this is an environment repair, not a changed physical experiment.
