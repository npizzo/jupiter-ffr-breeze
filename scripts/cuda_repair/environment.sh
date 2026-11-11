# Scoped CUDA preferences and compile cache; original package pins are unchanged.
export JULIA_LOAD_PATH="@:$FFR_ROOT/environments/cascade-cuda122:@stdlib"
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-cascade-cuda122:$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4 JULIA_NUM_PRECOMPILE_TASKS=4 OPENBLAS_NUM_THREADS=1
