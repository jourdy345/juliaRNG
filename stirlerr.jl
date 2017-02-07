function stirlerr(n::Float64)
    S0 = 0.083333333333333333333       
    S1 = 0.00277777777777777777778     
    S2 = 0.00079365079365079365079365  
    S3 = 0.000595238095238095238095238 
    S4 = 0.0008417508417508417508417508

    sferr_halves = [
        0.0, 
        0.1534264097200273452913848,  
        0.0810614667953272582196702,  
        0.0548141210519176538961390,  
        0.0413406959554092940938221,  
        0.03316287351993628748511048, 
        0.02767792568499833914878929, 
        0.02374616365629749597132920, 
        0.02079067210376509311152277, 
        0.01848845053267318523077934, 
        0.01664469118982119216319487, 
        0.01513497322191737887351255, 
        0.01387612882307074799874573, 
        0.01281046524292022692424986, 
        0.01189670994589177009505572, 
        0.01110455975820691732662991, 
        0.010411265261972096497478567, 
        0.009799416126158803298389475, 
        0.009255462182712732917728637, 
        0.008768700134139385462952823, 
        0.008330563433362871256469318, 
        0.007934114564314020547248100, 
        0.007573675487951840794972024, 
        0.007244554301320383179543912, 
        0.006942840107209529865664152, 
        0.006665247032707682442354394, 
        0.006408994188004207068439631, 
        0.006171712263039457647532867, 
        0.005951370112758847735624416, 
        0.005746216513010115682023589, 
        0.005554733551962801371038690
    ]
    M_LN_SQRT_2PI = 0.918938533204672741780329736406
    nn = 0.
    if n <= 15.
        nn = n+n
        if nn == floor(Int64,nn)
            return sferr_halves[floor(Int64,nn)]
        end
        return lgamma(n+1.)-(n+0.5)*log(n)+n-M_LN_SQRT_2PI
    end
    nn = n*n
    if n > 500
        return (S0-S1/nn)/n
    end
    if n > 80
        return (S0-(S1-S2/nn)/nn)/n
    end
    if n > 35
        return (S0-(S1-(S2-S3/nn)/nn)/nn)/n
    end
    return (S0-(S1-(S2-(S3-S4/nn)/nn)/nn)/nn)/n
end