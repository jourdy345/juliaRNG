function dgeom(k::Int64, p::Float64)
    if k < 0
        throw(ArgumentError("k should be greater than or equal to 0"))
    end
    if p <= 0. || p > 1.
        throw(ArgumentError("p should be 0 < p < 1"))
    end
    k_d = Float64(k)
    res = exp(k_d*log(1.-p)+log(p))
    res
end