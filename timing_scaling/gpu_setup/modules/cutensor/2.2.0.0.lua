-- -*- lua -*-

local name      = "cutensor"
local version   = "2.2.0.0"
whatis("Name         : " .. name)
whatis("Version      : " .. version)

family("cutensor")
--depends_on("cuda/12.1")

--conflict()

--In this section all dependencies will be specified

--binary location without architecture extension
local home    = "/u/area/ntosato/scratch/timing/time_scale_final/software/programs/cutensor/2.2.0.0"


prepend_path("CPATH", home .. "/include")
prepend_path("INCLUDE", home .. "/include")
prepend_path("LD_LIBRARY_PATH", home .. "/lib/12/")
prepend_path("LIBRARY_PATH", home .. "/lib/12/")
prepend_path("CUTENSOR_HOME", home)


