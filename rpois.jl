function rpois(n::Int64, λ::Float64)
    if n <= 0
        throw(ArgumentError("n should be greater than 0"))
    end
    if λ <= 0.
        throw(ArgumentError("λ should be greater than 0"))
    end
    res = zeros(n)
    for i = 1:n
        x = 0
        p = exp(-λ)
        s = p
        u = rand(1)[1]
        while u > s
            x += 1
            p *= λ/Float64(x)
            s += p
        end
        res[i] = x
    end
    res
end