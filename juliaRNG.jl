module juliaRNG
using Roots
export dexp,dgeom,dgamma,dnorm,pexp,pnorm,qnorm,rbeta,rchisq,rexp,rgamma,rgeom,rgig,rKS,rpois,rtexp,rtnorm,rWish
include("dexp.jl")
include("dgeom.jl")
include("dgamma.jl")
include("dnorm.jl")
include("pexp.jl")
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
include("stirlerr.jl")
include("bd0.jl")
end