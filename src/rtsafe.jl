export rtsafe


"""
    rtsafe(fdf, x0, a0, b0;
           [maxiter, xatol, xrtol, fatol])

Solve f(x) == 0.

The returned x is such that, if r is the actual root, then:

    abs(x - r) ≤ max(xatol, xrtol * abs(x))

or

    abs(f(x)) ≤ fatol

The function is specified by the argument fdf, which computes
the objective function and its derivative in a single call, as:

    f(x), f'(x) = fdf(x)

(a0,b0) bracket the root and must satisfy F(a0) * F(b0) ≤ 0.
x0 is an initial guess for the root and must satisfy

    x0 ∈ [a0, b0]
"""
function rtsafe(fdf,                    # f(x), f'(x) = fdf(x)
                x0, a0, b0;             # initial root guess and bracket
                maxiter::Integer = 100, # maximum number of iterations
                xatol::Real = 1e-10,    # absolute tolerance in abscisa
                xrtol::Real = 1e-10,    # relative tolerance in abscisa
                fatol::Real = 0,        # absolute tolerance in f(x)
                )

    #= Based on rtsafe, from Press, William H., et al. 2007. Numerical Recipes. 3rd edition.
    We use false position instead of bissection. At each iteration this routine does the best
    between a Newton step or a false position step. Thus it is guaranteed to converge to a root. =#

    @assert xatol ≥ 0 && xrtol ≥ 0 && fatol ≥ 0
    @assert maxiter ≥ 0
    @assert a0 ≤ x0 ≤ b0 || b0 ≤ x0 ≤ a0
    @assert isfinite(x0) && isfinite(a0) && isfinite(b0)

    a, b = minmax(a0, b0)   # current bracket
    x  = x0                 # current and last guess
    dx = abs(b - a)         # last step size and before last

    call_count = 0  # counts number of fdf calls
    F(x) = begin
        call_count += 1
        f, df = fdf(x)
        @assert !isnan(f) && !isnan(df)
        f, df
    end

    f, df = F(x)
    fa, dfa = F(a)
    fb, dfb = F(b)

    # bracket condition
    @assert fa ≤ 0 ≤ fb || fb ≤ 0 ≤ fa

    bissect_count = newton_count = 0

    for iter = 1 : maxiter
        xold  = x
        dxold = dx

        #= new root guess, either by Newton,
        or False Position =#

        newtondx = f / df
        newtonx = x - newtondx

        if a ≤ newtonx ≤ b && abs(2f) ≤ abs(dxold * df)
            # newton step good to go
            newton_count += 1
            x  = newtonx
            dx = abs(newtondx)
        else
            bissect_count += 1

            # do false position
            # x = a + (b-a)*fa/(fa-fb)
            # dx = max(b - x, x - a)

            # do bissection (seems faster than false position)
            #@info "bissection"
            dx = 0.5 * (b - a)
            x = a + dx

            if x == a
                @warn "convergence criteria cannot be met" x dx a b f df call_count bissect_count newton_count
                return x
            end
        end

        @assert a ≤ x ≤ b
        @assert dx ≥ 0

        f, df = F(x)

        # convergence criterion
        if abs(f) ≤ fatol || dx ≤ xatol || dx ≤ xrtol * abs(x)
            @debug "rtsafe converged" x xold dx a b f df call_count bissect_count newton_count
            @debug "rtsafe convergence criteria" xatol xrtol fatol (abs(f) ≤ fatol) (dx ≤ xatol) (dx ≤ xrtol * abs(x)) iter maxiter
            return x
        end

        # update bracket
        if sign(f) == sign(fa)
            a, fa, dfa = x, f, df
        elseif sign(f) == sign(fb)
            b, fb, dfb = x, f, df
        else
            error()
        end

        @assert a ≤ b
        @assert fa ≤ 0 ≤ fb || fb ≤ 0 ≤ fa
    end

    @error "Maximum number of iterations exceeded in rtsafe without convergence" x0 a0 b0 x a b f df dx fa fb call_count bissect_count newton_count
    @debug "rtsafe convergence criteria" xatol xrtol fatol (abs(f) ≤ fatol) (dx ≤ xatol) (dx ≤ xrtol * abs(x)) iter maxiter
    return x
end
