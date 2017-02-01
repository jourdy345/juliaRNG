function rexp(n::Int64, λ::Float64 = 1.)
    if n <= 0
        throw(ArgumentError("n should be greater than 0"))
    end
    if λ <= 0.
        throw(ArgumentError("λ should be greater than 0"))
    end
    res = -log(rand(n))/λ
    res
end