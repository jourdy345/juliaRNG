function rtexp(a::Float64, b::Float64)
    twoasq = 2.*a*a
    expab = expm1(-a*(b-a))
    z = 0.
    e = 0.
    while true
        z = log1p(rand(1)[1]*expab)
        e = -log(rand(1)[1])
        if twoasq*e <= z*z
            break
        end
    end
    return a-z/a
end