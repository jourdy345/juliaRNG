"""
DATE OF CREATION

    Feb.02.2017



AUTHOR
   
    Daeyoung Lim,
    Department of Statistics,
    Korea University



DESCRIPTION

    Generates truncated normal random variates.
    The algorithm used here is the accept-reject algorithm of Robert (1995)



REFERENCE

    Robert, C. P. (1995).
    Simulation of truncated normal variables. Statistics and computing,
    5(2), 121-125.
"""


function rtnorm(μ::Float64,σ::Float64,low::Float64,high::Float64)
    if low > high
        throw(ArgumentError("The lower bound should be smaller than the upper bound."))
    end
    draw = 0.
    typ = ""
    c_μ = μ
    c_σ = σ
    c_low = low
    c_high = high
    c_stdlow = (c_low-c_μ)/c_σ
    c_stdhigh = (c_high-c_μ)/c_σ

    if 0. <= c_stdhigh && 0. >= c_stdlow
        typ = "1"
    end
    if 0. < c_stdlow && c_stdhigh == Inf
        typ = "2"
    end
    if 0. > c_stdhigh && c_stdlow == -Inf
        typ = "3"
    end

    if (0. > c_stdhigh || 0. < c_stdlow) && !(isinf(c_stdhigh) || isinf(c_stdlow))
        typ = "4"
    end

    if typ == "1"
        z = 0.
        valid = false
        while !valid
            z = randn(1)[1]
            if c_stdlow <= z <= c_stdhigh
                valid = true
            end
        end
        draw = z
    end
    if typ == "3"
        c_stdlow = -1.*c_stdhigh
        c_stdhigh = Inf
        c_σ = -1.*c_σ
        typ = "2"
    end
    if typ == "2"
        αstar = (c_stdlow+sqrt(c_stdlow^2.+4.))/2.
        α = αstar
        e_ = 0.
        z = 0.
        ρ = 0.
        u = 0.
        valid = false
        while !valid
            e_ = rexp(1)[1]
            z = c_stdlow+e_/α
            ρ = exp((-α-z)^2./2.)
            u = rand(1)[1]
            if u <= ρ
                valid = true
            end
        end
        draw = z
    end
    if typ == "4"
        val1 = (2.*sqrt(exp(1.)))/(c_stdlow+sqrt(c_stdlow^2.+4.))
        val2 = val2 = exp((c_stdlow^2.-c_stdlow*sqrt(c_stdlow^2.+4.))/(4.))
        if c_stdhigh > c_stdlow+val1*val2
            valid = false
            while !valid
                αstar = (c_stdlow+sqrt(c_stdlow^2.+4.))/2.
                α = αstar
                e_ = 0.
                ρ = 0.
                u = 0.
                valid1 = false
                while !valid1
                    e_ = rexp(1)[1]
                    draw = c_stdlow+e_/α
                    ρ = exp((-α-draw)^2./2.)
                    u = rand(1)[1]
                    if u <= ρ
                        valid1 = true
                    end
                end
                if draw <= c_stdhigh
                    valid = true
                end
            end
        else
            valid = false
            ρ = 0.
            z = 0.
            u = 0.
            while !valid
                z = (c_stdhigh-c_stdlow)*rand(1)[1]+c_stdlow
                if 0 < c_stdlow
                    ρ = exp((c_stdlow^2.-z^2.)/2.)
                elseif c_stdhigh < 0.
                    ρ = exp((c_stdhigh^2.-z^2.)/2.)
                elseif 0. < c_stdhigh && c_stdlow < 0.
                    ρ = exp(-z^2./2.)
                end
                u = rand(1)[1]
                if u <= ρ
                    valid = true
                end
            end
            draw = z
        end
    end
    return c_μ+c_σ*draw
end