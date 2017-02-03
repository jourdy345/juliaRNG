"""
Date of Creation: Feb.02.2017
Author          : Daeyoung Lim,
                  Department of Statistics,
                  Korea University

Descriptiion    : This code generates random variables from
                  the generalized inverse Gaussian distribution
                  with the given parameters.
"""
function rgig(n::Int64, λ::Float64, χ::Float64, ψ::Float64)
    if χ < 0. || ψ < 0.
        throw(ArgumentError("Invalid parameters for GIG"))
    end
    if χ == 0. && λ <= 0.
        throw(ArgumentError("Invalid parameters for GIG")) 
    end
    if ψ == 0. && λ >= 0.
        throw(ArgumentError("Invalid parameters for GIG"))
    end
    if χ == 0. && λ > 0.
        return rgamma(n,λ,2./ψ)
    end
    if ψ == 0. && λ < 0.
        return 1./rgamma(n,-λ,2./χ)
    end

    if λ == -0.5
        U = rand(n)
        V = randn(n).^2.
        δ = sqrt(χ)
        γ = sqrt(ψ)

        z1 = δ/γ+V./(2.*γ^2.)-sqrt(V.*δ/(γ^3.)+(V./(2.*γ^2.)).^2.)
        z2 = ((δ/γ)^2.)./z1
        pz1 = δ./(δ+γ*z1)
        s = (1.-sign(U-pz1))./2.
        return z1.*s+z2.*(1.-s)
    elseif λ == 1.
        α = sqrt(ψ/χ)
        β = sqrt(ψ*χ)
        m = sign(β)
        g1(y) = 0.5*β*y^3.-y^2.*(0.5*β+2.)+y*(-0.5*β)+0.5*β
        upper = 1.
        while g1(upper) <= 0.
            upper *= 2.
        end
        yM = fzero(g1,0.,1.)
        yP = fzero(g1,1.,upper)
        a = (yP-1.)*exp(-0.25*β*(yP+1./yP-2.))
        b = (yM-1.)*exp(-0.25*β*(yM+1./yM-2.))
        c = -0.25*β*2.
        output = zeros(n)
        for i = 1:n
            Y = 0.
            needValue = true
            while needValue
                R1 = rand(1)[1]
                R2 = rand(1)[1]
                Y = 1.+a*R2/R1+b*(1.-R2)/R1
                if Y > 0.
                    if -log(R1) >= 0.25*β*(Y+1./Y)+c
                        needValue = false
                    end
                end
            end
            output[i] = Y
        end
        return output/α
    else
        α = sqrt(ψ/χ)
        β = sqrt(ψ*χ)
        m = ((λ-1.)+sqrt((λ-1.)^2.+β^2.))/β
        g2(y) = 0.5*β*y^3.-y^2.*(0.5*β*m+λ+1.)+y*((λ-1.)*m-0.5*β)+0.5*β*m
        upper = m
        while g2(upper) <= 0.
            upper *= 2.
        end
        yM = fzero(g2,0,m)
        yP = fzero(g2,m,upper)
        a = (yP-m)*(yP/m)^(0.5*(λ-1.))*exp(-0.25*β*(yP+1./yP-m-1./m))
        b = (yM-m)*(yM/m)^(0.5*(λ-1.))*exp(-0.25*β*(yM+1./yM-m-1./m))
        c = -0.25*β*(m+1./m)+0.5*(λ-1.)*log(m)
        output = zeros(n)
        for i = 1:n
            needValue = true
            Y = 0.
            while needValue
                R1 = rand(1)[1]
                R2 = rand(1)[1]
                Y = m+a*R2/R1+b*(1.-R2)/R1
                if Y > 0.
                    if -log(R1) >= -0.5*(λ-1.)*log(Y)+0.25*β*(Y+1./Y)+c
                        needValue = false
                    end
                end
            end
            output[i] = Y
        end
        return output/α
    end
end









