include("rchisq.jl")
function rWish(df::Int64, S::Array{Float64,2})
    if size(S,1) != size(S,2)
        throw(ArgumentError("S should be square symmetric positive definite"))
    end
    m = size(S,1)
    if df <= m-1
        throw(ArgumentError("The degree(s) of freedom should be greater than nrow(S)-1"))
    end
    A = zeros(m,m)
    for i = 1:m
        A[i,i] = sqrt(rchisq(1,Float64(df-i)+1.)[1])
    end
    for j = 1:(m-1)
        for i = (j+1):m
            A[i,j] = randn(1)[1]
        end
    end
    C = At_mul_B(A,chol(S))
    C = At_mul_B(C,C)
    C
end