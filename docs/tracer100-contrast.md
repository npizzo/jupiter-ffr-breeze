# Tracer contrast for the completed 100-day movies

Requested November 11, 2026: make late-time passive-tracer structures visible in the full 0–100-day heated/unheated movies, and publish the revised movies to the GitHub gallery. This changes only rendering. Original movies and all simulation fields are retained.

Unity audit job 1437416 found that the tracer is mixed toward a background concentration near C = 0.0793. At day 100 the heated lower/middle slices still contain median absolute deviations of about 0.0018 and 90th-percentile deviations of about 0.0062. The Float32 spacing near the mean is about 7.45e-9; these variations are well above saved precision. The original broad zero-centered sequential scale gave those physically distinct concentrations nearly identical colors.

The contrast variant uses the diverging brown–teal `BrBG` colormap with a fixed center C0 = 0.0793082 and a near-center asinh width 0.001. The negative and positive transform branches are normalized to map the common background to the colormap midpoint while retaining the full observed concentration range (including 0–1). Physical colorbar ticks remain actual concentration C. Data are not mean-subtracted, clipped, or rescaled per frame. The same mapping is used across both runs, all three depths, and all 801 saved times.

Vorticity panels retain their original fixed scales; rendering checks them against the original full-100-day metadata. The three requested outputs are the matched comparison, heated case, and unheated case. Rendering array 1437852 runs on Unity, followed by complete movie decoding and frame/time checks within each task. Only finished publication assets are transferred for GitHub upload.

This is a visualization improvement, not evidence of numerical convergence or a new transport result. The original absolute-concentration videos remain available for comparison. [Watch on the GitHub gallery](https://npizzo.github.io/jupiter-ffr-breeze/#tracer100-contrast).
