@inline cheb_scale(x, bounds) = 2(x-bounds[1])/(bounds[2]-bounds[1]) - 1
@inline cheb_unscale(x, bounds) = (x+1)*(bounds[2]-bounds[1])/2 + bounds[1]


function apply_jackson_kernel!(coefs)
    Mp = lastindex(coefs) + 2
    for i in eachindex(coefs)
        m = i-1
        coefs[i] *= (1/Mp)*((Mp-m)cos(m*π/Mp) + sin(m*π/Mp)/tan(π/Mp))
    end
    return coefs
end


function cheb_coefs!(M, func, bounds; buf, plan)
    N = length(buf)
    for i in eachindex(buf)
        x_i = cos((i-0.5)π / N)
        buf[i] = func(cheb_unscale(x_i, bounds))
    end

    mul!(buf, plan, buf)
    buf ./= N
    buf[1] /= 2
    return view(buf, 1:M)
end

"""
    cheb_coefs(M, nsamples, func, bounds)

Generate `M` coefficients of the Chebyshev expansion using `nsamples` of the
function `func`. Sample points are taken within the interval specified by
`bounds = (lo, hi)`.
"""
function cheb_coefs(M, nsamples, func, bounds)
    @assert nsamples >= M
    buf = zeros(nsamples)
    plan = FFTW.plan_r2r!(buf, FFTW.REDFT10)
    return cheb_coefs!(M, func, bounds; buf, plan)
end

"""
    cheb_eval(x, bounds, coefs)

Evaluate a function, specified in terms of the Chebyshev `coefs`, at point `x`.
`bounds` specifies the domain of the function.
"""
function cheb_eval(x, bounds, coefs)
    x = cheb_scale(x, bounds)
    Tn2, Tn1 = 1, x
    ret = coefs[1]*Tn2 + coefs[2]*Tn1
    for coef in @view coefs[3:end]
        Tn3 = 2x*Tn1 - Tn2
        ret += coef*Tn3
        Tn1, Tn2 = Tn3, Tn1
    end
    return ret
end
