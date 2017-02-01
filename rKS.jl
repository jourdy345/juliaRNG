function rKS(n::Int64)
    t = 0.75
    tprime = pi^2./(8.*t^2.)
    p = 0.373
    q = 1.-p
    output = zeros(n)
    for i = 1:n
        X = 0.
        if rand(1)[1] < exp(log(p)-log(p+q))
            while true
                accept = false
                G = 0.
                while !accept
                    E = rexp(2)
                    E[1] /= 1.-0.5/tprime
                    E[2] *= 2.
                    G = tprime+E[1]
                    accept = (E[1]^2. <= tprime*E[2]*(G+tprime))
                    if !accept
                        accept = (G/tprime-1.-log(G/tprime) <= E[2])
                    end
                end

                ## Devroye is clever
                X = pi/sqrt(8.*G)
                W = 0.
                Z = 0.5/G
                P = exp(-G)
                n = 1.
                Q = 1.
                U = rand(1)[1]
                go = true
                success = 0
                while go
                    W += Z*Q
                    if U >= W
                        output[i] = X
                        success = 1
                        break
                    end
                    n += 2.
                    Q = P^(n^2.-1.)
                    W -= n^2.*Q
                    go = U >= W
                end
                if success == 1
                    break
                end
            end
        else
            while true
                E = rexp(1)[1]
                U = rand(1)[1]
                X = sqrt(t^2.+0.5*E)
                W = 0.
                n = 1.
                Z = exp(-2.*X^2.)
                go = true
                success = 0
                while go
                    n += 1.
                    W += n^2.*Z^(n^2.-1.)
                    if U >= W
                        output[i] = X
                        success = 1
                        break
                    end
                    n += 1.
                    W -= n^2.*Z^(n^2.-1)
                    go = U >= W
                end
                if success == 1
                    break
                end
            end
        end
    end
    output
end

n = 0
while true
    n += 1
    if n % 2 == 0
        continue
    end
    println(n)
    if n > 15
        break
    end
end