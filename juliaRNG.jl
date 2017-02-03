module juliaRNG
using Roots
export dgeom,dnorm,pnorm,qnorm,rbeta,rchisq,rexp,rgamma,rgeom,rgig,rKS,rpois,rtexp,rtnorm,rWish
include("dgeom.jl")
include("dnorm.jl")
include("pnorm.jl")
include("qnorm.jl")
include("rbeta.jl")
include("rchisq.jl")
include("rexp.jl")
include("rgamma.jl")
include("rgeom.jl")
include("rgig.jl")
include("rKS.jl")
include("rpois.jl")
include("rtexp.jl")
include("rtnorm.jl")
include("rWish.jl")
end