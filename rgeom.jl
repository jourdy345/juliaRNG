function rgeom(n::Int64, p::Float64)
    if n <= 0
        throw(ArgumentError("n should be greater than 0"))
    end
    if p <= 0. || p > 1.
        throw(ArgumentError("p should be 0 < p < 1"))
    end
    res = zeros(n)
    for i = 1:n
        X = 0
        S = dgeo(0,p)
        u = rand(1)[1]
        while u > S
            X += 1
            S += dgeom(X,p)
        end
        res[i] = X
    end
    res
end