"""
DATE OF CREATION
    
    Feb.02.2017


AUTHOR
   
    Daeyoung Lim,
    Department of Statistics,
    Korea University


REFERENCE

    R. C. H. Cheng (1978).
    Generating beta variates with nonintegral shape parameters.
    Communications of the ACM 21, 317-322.
    (Algorithms BB and BC)

"""
function rbeta(aa::Float64,bb::Float64)
    if aa < 0. || bb < 0.
        return NaN
    end
    if isinf(aa) && isinf(bb)
        return 0.5
    end
    if aa == 0. && bb == 0.
        return rand(1)[1] < 0.5 ? 0. : 1.
    end
    if isinf(aa) || bb == 0.
        return 1.
    end
    if isinf(bb) || aa == 0.
        return 0.
    end

    DBL_MAX = 1.7976931348623157e+308
    DBL_MAX_EXP = 1024
    M_LN2 = 0.69314718055994530941723212146
    expmax = DBL_MAX_EXP*M_LN2
    olda = -1.
    oldb = -1.
    β = 0.
    δ = 0.
    γ = 0.
    k1 = 0.
    k2 = 0.
    r = 0.
    s = 0.
    t = 0.
    u1 = 0.
    u2 = 0.
    v = 0.
    w = 0.
    y = 0.
    z = 0.
    qsame = (olda == aa) && (oldb == bb)
    if !qsame
        olda = aa
        oldb = bb
    end

    a = (isnan(aa) || isnan(bb)) ? aa+bb : ((aa < bb) ? aa : bb)
    b = (isnan(aa) || isnan(bb)) ? aa+bb : ((aa < bb) ? bb : aa)

    α = a+b
    if a <= 1.
        if !qsame
            β = 1./a
            δ = 1.+b-a
            k1 = δ*(0.0138889+0.0416667*a)/(b*β-0.777778)
            k2 = 0.25+(0.5+0.25/δ)*a
        end
        while true
            u1 = rand(1)[1]
            u2 = rand(1)[1]
            if u1 < 0.5
                y = u1*u2
                z = u1*y
                if 0.25*u2+z-y >= k1
                    continue
                end
            else
                z = u1*u1*u2
                if z <= 0.25
                    v = β*log(u1/(1.-u1))
                    if v <= expmax
                        w = b*exp(v)
                        if isinf(w)
                            w = DBL_MAX
                        end
                    else
                        w = DBL_MAX
                    end
                    break
                end
                if z >= k2
                    continue
                end
            end
            v = β*log(u1/(1.-u1))
            if v <= expmax
                w = b*exp(v)
                if isinf(w)
                    w = DBL_MAX
                end
            else
                w = DBL_MAX
            end

            if α*(log(α/(a+w))+v)-1.3862944 >= log(z)
                break
            end
        end
        return aa == a ? a/(a+w) : w/(a+w)
    else
        if !qsame
            β = sqrt((α-2.)/(2.*a*b-α))
            γ = a+1./β
        end
        while true
            u1 = rand(1)[1]
            u2 = rand(1)[1]
            v = β*log(u1/(1.-u1))
            if v <= expmax
                w = a*exp(v)
                if isinf(w)
                    w = DBL_MAX
                end
            else
                w = DBL_MAX
            end
            z = u1*u1*u2
            r = γ*v-1.3862944
            s = a+r-w
            if s+2.609438 >= 5.0*z
                break
            end
            t = log(z)
            if s > t
                break
            end
            if r+α*log(α/(b+w)) >= t
                break
            end
        end
        return aa != a ? b/(b+w) : w/(b+w)
    end
end