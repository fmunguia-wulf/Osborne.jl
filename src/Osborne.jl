module Osborne

import IMASdd
import IMAS
import EFIT
using FastInterpolations

include("io.jl")
export readp, PFile, PFileProfile, PFileIonSpecies, interpolate_1d_profile

include("pfile_imas.jl")
export pfile2imas!

end
