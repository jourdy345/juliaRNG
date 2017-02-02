"""
DATE OF CREATION
    
    Feb.02.2017



AUTHOR  
  
    Daeyoung Lim,
    Department of Statistics,
    Korea University



DESCRIPTION
 
    Compute the quantile function for the normal distribution.
    For small to moderate probabilities, algorithm referenced
    below is used to obtain an initial approximation which is
    polished with a final Newton step.
  
    For very large arguments, an algorithm of Wichura is used.
 

 
REFERENCE
 
    Beasley, J. D. and S. G. Springer (1977).
    Algorithm AS 111: The percentage points of the normal distribution,
    Applied Statistics, 26, 118-121.
  
    Wichura, M.J. (1988).
    Algorithm AS 241: The Percentage Points of the Normal Distribution.
    Applied Statistics, 37, 477-484.
"""

function qnorm(p::Float64, μ::Float64, σ::Float64, lower_tail::Bool, log_p::Bool)
    if isnan(p) || isnan(μ) || isnan(σ)
        return NaN
    end
    # R_D_0 = log_p ? -Inf : 0.
    # R_D_1 = log_p ? 0. : 1.
    # R_DT_0 = lower_tail ? R_D_0 : R_D_1
    # R_DT_1 = lower_tail ? R_D_1 : R_D_0

    if log_p
        if p > 0.
            return NaN
        end
        if p == 0.
            return lower_tail ? Inf : -Inf
        end
        if p == -Inf
            return lower_tail ? -Inf : Inf
        end
    else
        if p < 0. || p > 1.
            return NaN
        end
        if p == 0.
            return lower_tail ? -Inf : Inf
        end
        if p == 1.
            return lower_tail ? Inf : -Inf
        end
    end

    if σ < 0.
        return NaN
    end
    if σ == 0.
        return μ
    end
    p_ = log_p ? (lower_tail ? exp(p) : -expm1(p)) : (lower_tail ? (p) : (0.5-(p)+0.5))
    q = p_-0.5
    r = 0.
    val = 0.
    if abs(q) <= .425
        r = .180625-q*q
        val = q*(((((((r*2509.0809287301226727+33430.575583588128105)*r+67265.770927008700853)*r+45921.953931549871457)*r+13731.693765509461125)*r+1971.5909503065514427)*r+133.14166789178437745)*r+3.387132872796366608)/(((((((r*5226.495278852854561+28729.085735721942674)*r+39307.89580009271061)*r+21213.794301586595867)*r+5394.1960214247511077)*r+687.1870074920579083)*r+42.313330701600911252)*r+1.)
    else
        if q > 0.
            r = log_p ? (lower_tail ? -expm1(p) : exp(p)) : (lower_tail ? (0.5-(p)+0.5) : (p))
        else
            r = p_
        end

        r = sqrt(- ((log_p && ((lower_tail && q <= 0.) || (!lower_tail && q > 0.))) ? p : log(r)))

        if r <= 5.
            r += -1.6
            val = (((((((r*7.7454501427834140764e-4+.0227238449892691845833)*r+.24178072517745061177)*r+1.27045825245236838258)*r+3.64784832476320460504)*r+5.7694972214606914055)*r+4.6303378461565452959)*r+1.42343711074968357734)/(((((((r*1.05075007164441684324e-9+5.475938084995344946e-4)*r+.0151986665636164571966)*r+.14810397642748007459)*r+.68976733498510000455)*r+1.6763848301838038494)*r+2.05319162663775882187)*r+1.)
        else
            r += -5.
            val = (((((((r*2.01033439929228813265e-7+2.71155556874348757815e-5)*r+.0012426609473880784386)*r+.026532189526576123093)*r+.29656057182850489123)*r+1.7848265399172913358)*r+5.4637849111641143699)*r+6.6579046435011037772)/(((((((r*2.04426310338993978564e-15+1.4215117583164458887e-7)*r+1.8463183175100546818e-5)*r+7.868691311456132591e-4)*r+.0148753612908506148525)*r+.13692988092273580531)*r+.59983220655588793769)*r+1.)
        end

        if q < 0.
            val *= -1
        end
    end
    return μ+σ*val
end