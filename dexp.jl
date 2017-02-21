function dexp(x::Float64,scale::Float64,give_log::Bool)
    if isnan(x) || isnan(scale)
        return x+scale
    end
    if scale <= 0.
        return NaN
    end
    R_D_0 = give_log ? -Inf : 0.
    R_D_1 = give_log ? 0. : 1.
    if x < 0.
        return R_D_0
    end
    return give_log ? (-x/scale)-log(scale) : exp(-x/scale)/scale
end