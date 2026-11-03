module MicroscopePSFsForwardDiffExt

import ForwardDiff
import MicroscopePSFs: _primal

_primal(d::ForwardDiff.Dual) = _primal(ForwardDiff.value(d))

end
