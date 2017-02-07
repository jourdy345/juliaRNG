function dgamma(x::Float64,shape::Float64,scale::Float64,give_log::Bool)
    if isnan(x) || isnan(shape) || isnan(scale)
        return x+shape+scale
    end

    if shape < 0. || scale <= 0.
        return NaN
    end

    R_D_0 = give_log ? -Inf : 0.
    R_D_1 = give_log ? 0. : 1.
    DBL_MIN = 2.22507385850720140000e-308
    M_2PI = 6.2831853071795865
    if x < 0.
        return R_D_0
    end
    if shape == 0.
        return x == 0 ? Inf : R_D_0
    end
    if x == 0.
        if shape < 1.
            return Inf
        end
        if shape > 1.
            return R_D_0
        end
        return give_log ? -log(scale) : 1./scale
    end

    pr = 0.
    if shape < 1.
        if x/scale == 0.
            pr = shape == 0. ? R_D_1 : R_D_0
        end
        if isinf(x/scale)
            pr = R_D_0
        end
        if shape < 0.
            pr = R_D_0
        end
        if shape <= x/scale*DBL_MIN
            pr = give_log ? -x/scale+shape*log(x/scale)-lgamma(shape+1.) : exp(-x/scale+shape*log(x/scale)-lgamma(shape+1.))
        end
        if x/scale < shape*DBL_MIN
            pr = give_log ? -x/scale+shape*log(x/scale)-lgamma(x+1.) : exp(-x/scale+shape*log(x/scale)-lgamma(x+1.))
        end
        
        pr = give_log ? -0.5*log(M_2PI*shape)-stirlerr(shape)-bd0(shape,x/scale) : exp(-stirlerr(shape)-bd0(shape,x/scale))/sqrt(M_2PI*shape)
        return give_log ? pr+log(shape/x) : pr*shape/x
    end

    if x/scale == 0.
        pr = shape-1. == 0. ? R_D_1 : R_D_0
    end
    if isinf(x/scale)
        pr = R_D_0
    end
    if shape-1. < 0.
        pr = R_D_0
    end
    if shape-1. <= x/scale*DBL_MIN
        pr = give_log ? -x/scale+(shape-1.)*log(x/scale)-lgamma((shape-1.)+1.) : exp(-x/scale+(shape-1.)*log(x/scale)-lgamma((shape-1.)+1.))
    end
    if x/scale < (shape-1.)*DBL_MIN
        pr = give_log ? -x/scale+(shape-1.)*log(x/scale)-lgamma(x+1.) : exp(-x/scale+(shape-1.)*log(x/scale)-lgamma(x+1.))
    end
    pr = give_log ? -0.5*log(M_2PI*(shape-1.))-stirlerr((shape-1.))-bd0((shape-1.),x/scale) : exp(-stirlerr((shape-1.))-bd0((shape-1.),x/scale))/sqrt(M_2PI*(shape-1.))
    return give_log ? pr - log(scale) : pr/scale
end

"""
double attribute_hidden dpois_raw(double x, double lambda, int give_log)
{
    if (lambda == 0) return( (x == 0) ? R_D__1 : R_D__0 );
    if (!R_FINITE(lambda)) return R_D__0;
    if (x < 0) return( R_D__0 );
    if (x <= lambda * DBL_MIN) return(R_D_exp(-lambda) );
    if (lambda < x * DBL_MIN) return(R_D_exp(-lambda + x*log(lambda) -lgammafn(x+1)));
    return(R_D_fexp( M_2PI*x, -stirlerr(x)-bd0(x,lambda) ));
}
"""









