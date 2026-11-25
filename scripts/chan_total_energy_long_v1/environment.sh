# Private package environment; active experiments retain their original pins.
export TE_PROJECT="$FFR_ROOT/environments/chan-total-energy-8273be8"
export JULIA_LOAD_PATH="@:@stdlib"
export JULIA_DEPOT_PATH="$FFR_ROOT/depot-chan-total-energy-8273be8:$FFR_ROOT/depot-cascade-cuda122:$FFR_ROOT/depot"
export JULIA_NUM_THREADS=4 JULIA_NUM_PRECOMPILE_TASKS=4 OPENBLAS_NUM_THREADS=1
export JULIA_PKG_PRECOMPILE_AUTO=0
