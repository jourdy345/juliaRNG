"""
DATE OF CREATION
    
    Feb.02.2017


AUTHOR
   
    Daeyoung Lim,
    Department of Statistics,
    Korea University

DESCRIPTION

    Compute the density of the normal distribution
"""
function dnorm(x::Float64,μ::Float64,σ::Float64,give_log::Bool)
    if isnan(x) || isnan(μ) || isnan(σ)
        return NaN
    end
    R_D_0         = give_log ? -Inf : 0.
    R_D_1         = give_log ? 0. : 1.
    DBL_MAX       = 1.7976931348623157e+308
    DBL_MIN_EXP   = -1022
    DBL_MANT_DIG  = 53
    M_LN_SQRT_2PI = 0.918938533204672741780329736406
    M_LN2         = 0.69314718055994530941723212146
    M_1_SQRT_2PI  = 0.398942280401432677939946059934

    if isinf(σ)
        return R_D_0
    end
    if isinf(x) && μ == x
        return NaN
    end
    if σ <= 0.
        if σ < 0.
            throw(ArgumentError("σ should be non-negative"))
        end
        return x == μ ? Inf : R_D_0
    end
    x = (x-μ)/σ
    if isinf(x)
        return R_D_0
    end
    x = abs(x)
    if x >= 2.*sqrt(DBL_MAX)
        return R_D_0
    end
    if give_log
        return -(M_LN_SQRT_2PI+0.5*x*x+log(σ))
    end

    if x < 5.
        return M_1_SQRT_2PI*exp(-0.5*x*x)/σ
    end
    if x > sqrt(-2.*M_LN2*(DBL_MIN_EXP+1.-DBL_MANT_DIG))
        return 0.
    end
    x1 = ldexp(Int64(ldexp(x,16)),-16)
    x2 = x-x1
    return M_1_SQRT_2PI/σ*(exp(-0.5*x1*x1)*exp((-0.5*x2-x1)*x2))
end