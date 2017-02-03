function rchisq(n::Int64, df::Float64)
    if isinf(df) || df < 0.0
        throw(ArgumentError("invalid value for df"))
    end
    output = rgamma(n,df/2.,2.)
    output
end