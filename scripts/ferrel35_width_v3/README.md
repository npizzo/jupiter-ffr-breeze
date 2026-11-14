# 35N native-budget restart and width sensitivity

Isolated copy of the existing 35N model. Physics is unchanged in the diagnostic day200–210 branch. The second experiment doubles north–south extent and cell count, preserving grid spacing, heating, cap treatment, f0/beta, and initial eddy speed/wavelength statistics. Initial realizations differ because the channel eigenmodes differ; the old checkpoint is not stretched.

Every complete SSPRK3 step records and integrates the model's native zonal momentum terms. Detailed step files are sparse; cumulative integrals accompany every three-hour movie frame. The centered mean/eddy decomposition and its difference from native advection remain explicit. This difference is not declared to be physical eddy driving. No changes are made to active campaign source or the CUDA environment.

Preflight requires unchanged physical parameters, continuous/restarted agreement (including integrated momentum terms), a mature day200 comparison against the original driver, and a full-size 256×512×64 six-hour GPU probe. The wider run is also gated on successful completion and numerical budget review of the ten-day diagnostic branch. Twenty sequential ten-day jobs reach day200 in the wide channel. All computations, analysis, and rendering belong on Unity.

Saved outputs retain actual vorticity and tracer at three depths, zonal heat/tracer fluxes, mass circulation, and native plus approximate mean/eddy momentum diagnostics. The wide-channel budget spans each ten-day segment; add integrated arrays across segment endpoints for longer intervals. Diagnostic sums intentionally restart from zero and do not modify physical state.
