-- -*- lua -*-

local name      = "R"
local version   = "4.3.3-h100"
whatis("Name         : " .. name)
whatis("Version      : " .. version)

family("R")
--depends_on("openBLAS")

--conflict()

--In this section all dependencies will be specified

--binary location without architecture extension
local home    = "/u/area/ntosato/scratch/timing/time_scale_final/software/programs/R/4.3.3-h100"


prepend_path("PATH", home .. "/bin")
prepend_path("LD_LIBRARY_PATH", home .."/lib64")
prepend_path("MANPATH", home .."/share/man")
prepend_path("R_LIBS_USER","/u/area/ntosato/scratch/timing/time_scale_final/software/programs/r_packages_h100")

