"""
DATE OF CREATION
    
    Feb.02.2017


AUTHOR
   
    Daeyoung Lim,
    Department of Statistics,
    Korea University


DESCRIPTION

    Pseudo-random variates from a t distribution


NOTES

    This function calls rchisq and randn to do the real work
"""
function rt(n::Int64, df::Float64)
    if isnan(df) || df <= 0.
        return NaN
    end

    output = zeros(n)
    for i = 1:n
        if isinf(df)
            output[i] = randn(1)[1]
        else
            num = randn(1)[1]
            output[i] = num/sqrt(rchisq(1,df)[1]/df)
        end
    end
    return output
end