"""
DATE OF CREATION
    
    Feb.02.2017



AUTHOR  
  
    Daeyoung Lim,
    Department of Statistics,
    Korea University


DESCRIPTION

    Random variates from the gamma distribution.



REFERENCES

    [1] Shape parameter a >= 1.  Algorithm GD in:

    Ahrens, J.H. and Dieter, U. (1982).
    Generating gamma variates by a modified
    rejection technique.
    Comm. ACM, 25, 47-54.


    [2] Shape parameter 0 < a < 1. Algorithm GS in:

    Ahrens, J.H. and Dieter, U. (1974).
    Computer methods for sampling from gamma, beta,
    poisson and binomial distributions.
    Computing, 12, 223-246.
"""

function rgamma(n::Int64, a::Float64, scale::Float64)
    if isinf(a) || isinf(scale) || a<0. || scale <= 0.
        throw(ArgumentError("Invalid parameter(s)"))
    end
    sqrt32 = 5.656854
    exp_m1 = 0.36787944117144232159
    q1 = 0.04166669
    q2 = 0.02083148
    q3 = 0.00801191
    q4 = 0.00144121
    q5 = -7.388e-5
    q6 = 2.4511e-4
    q7 = 2.424e-4

    a1 = 0.3333333
    a2 = -0.250003
    a3 = 0.2000062
    a4 = -0.1662921
    a5 = 0.1423657
    a6 = -0.1367177
    a7 = 0.1233795


    res = zeros(n)
    for i = 1:n
        aa = 0.
        aaa = 0.
        s = 0.
        s2 = 0.
        d = 0.
        q0 = 0.
        b = 0.
        si = 0.
        c = 0.
        e_ = 0.
        p = 0.
        q = 0.
        r = 0.
        t = 0.
        u = 0.
        v = 0.
        w = 0.
        x = 0.
        ret_val = 0.
        if a < 1.
            if a == 0.
                res[i] = 0.
                continue
            end
            e_ = 1.+exp_m1*a
            while true
                p = e_*rand(1)[1]
                if p >= 1.
                    x = -log((e_-p)/a)
                    if rexp(1,1.)[1] >= (1.-a)*log(x)
                        break
                    end
                else
                    x = exp(log(p)/a)
                    if rexp(1,1.)[1] >= x
                        break
                    end
                end
            end
            res[i] = scale*x
            continue
        end

        if a != aa
            aa = a
            s2 = a-0.5
            s = sqrt(s2)
            d = sqrt32-s*12.
        end
        t = randn(1)[1]
        x = s+0.5*t
        ret_val = x*x
        if t >= 0.
            res[i] = scale*ret_val
            continue
        end
        u = rand(1)[1]
        if d*u <= t*t*t
            res[i] = scale*ret_val
            continue
        end

        if a != aaa
            aaa = a
            r = 1./a
            q0 = ((((((q7*r+q6)*r+q5)*r+q4)*r+q3)*r+q2)*r+q1)*r

            if a <= 3.686
                b = 0.463+s+0.178*s2
                si = 1.235
                c = 0.195/s-0.079+0.16*s
            elseif a <= 13.022
                b = 1.654+0.0076*s2
                si = 1.68/s+0.275
                c = 0.062/s+0.024
            else
                b = 1.77
                si = 0.75
                c = 0.1515/s
            end
        end
        if x > 0.
            v = t/(s+s)
            if abs(v) <= 0.25
                q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
            else
                q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v)
            end
            if log(1.0-u) <= q
                res[i] = scale*ret_val
                continue
            end
        end
        while true
            e_ = rexp(1,1.)[1]
            u = rand(1)[1]
            u = u+u-1.
            if u < 0.
                t = b-si*e_
            else
                t = b+si*e_
            end
            if t >= -0.71874483771719
                v = t/(s+s)
                if abs(v) <= 0.25
                    q = q0+0.5*t*t*((((((a7*v+a6)*v+a5)*v+a4)*v+a3)*v+a2)*v+a1)*v
                else
                    q = q0-s*t+0.25*t*t+(s2+s2)*log(1.0+v)
                end
                if q > 0.
                    w = expm1(q)
                    if c*abs(u) <= w*exp(e_ - 0.5*t*t)
                        break
                    end
                end
            end
        end
        x = s+0.5*t
        res[i] = scale*x*x
    end
    res
end