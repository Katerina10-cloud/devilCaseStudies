-- -*- lua -*-

local name      = "openBLAS"
local version   = "0.3.29-a100"

whatis("Name         : " .. name)
whatis("Version      : " .. version)

family("BLAS")

local home    = "/u/area/ntosato/scratch/timing/time_scale_final/software/programs/openBLAS/0.3.29-a100/"

prepend_path{"PATH", home .. "/bin",delim=":",priority="0"}
prepend_path{"LD_LIBRARY_PATH", home .. "/lib",delim=":",priority="0"}
prepend_path{"LIBRARY_PATH", home .. "/lib",delim=":",priority="0"}
prepend_path{"CPATH", home .. "/include",delim=":",priority="0"}
setenv("OPENBLAS_DIR", home )
setenv("OPENBLAS_ROOT", home)
setenv("OPENBLAS_LIB", home .. "/lib")
setenv("OPENBLAS_IN",home .. "/include")
