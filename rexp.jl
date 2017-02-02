
"""
DATE OF CREATION
    
    Feb.02.2017



AUTHOR  
  
    Daeyoung Lim,
    Department of Statistics,
    Korea University



DESCRIPTION
 
    Random variates from the exponential distribution.
    Uses probability integral transform from the inverse CDF
    of the exponential distribution.
"""
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