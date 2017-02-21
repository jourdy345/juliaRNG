function pexp(x::Float64,scale::Float64,lower_tail::Bool,log_p::Bool)
    if isnan(x) || isnan(scale)
        return x+scale
    end
    if scale <= 0.
        return NaN
    end
    R_D_0 = log_p ? -Inf : 0.
    R_D_1 = log_p ? 0. : 1.
    R_DT_0 = lower_tail ? R_D_0 : R_D_1
    R_DT_1 = lower_tail ? R_D_1 : R_D_0
    if x <= 0.
        return R_DT_0
    end
    x = -(x/scale)
    M_LN2 = 0.693147180559945309417232121458
    return lower_tail ? (log_p ? (x > -M_LN2 ? log(-expm1(x)) : log1p(-exp(x))) : -expm1(x)) : (log_p ? x : exp(x))
end