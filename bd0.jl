function bd0(x::Float64,np::Float64)
    if isinf(x) || isinf(np) || np == 0.
        return NaN
    end
    ej = 0.
    s = 0.
    s1 = 0.
    v = 0.
    if abs(x-np) < 0.1*(x+np)
        v = (x-np)/(x+np)
        s = (x-np)*v
        ej = 2.*x*v
        v = v*v
        j = 1
        while true
            ej *= v
            s1 = s+ej/((j<<1)+1)
            if s1 == s
                return s1
            end
            s = s1
        end
    end
    return x*log(x/np)+np-x
end