function testrtsafe(fdf, eqn; x0, xl, xh, root)
    @info "solving " * eqn * " from x0 = $x0 (real root = $root) ... "
    x = rtsafe(fdf, x0, xl, xh)
    @test x ≈ root  atol=1e-15
end


@testset "rtsafe" begin
    fdf(x) = (exp(x) - ℯ, exp(x))
    testrtsafe(fdf, "exp(x) = e"; x0=5, xl=-10, xh=10, root=1)

    fdf(x) = (exp(x) - ℯ, exp(x))
    testrtsafe(fdf, "exp(x) = e"; x0=5, xl=10, xh=-10, root=1)

    fdf(x) = (exp(-x) - ℯ, -exp(-x))
    testrtsafe(fdf, "exp(-x) = e"; x0=5, xl=-10, xh=10, root=-1)

    fdf(x) = (exp(-x) - ℯ, -exp(-x))
    testrtsafe(fdf, "exp(-x) = e"; x0=5, xl=10, xh=-10, root=-1)


    fdf(x) = (exp(x) - ℯ, exp(x))
    testrtsafe(fdf, "exp(x) = e"; x0=-5, xl=-10, xh=10, root=1)

    fdf(x) = (sin(x), cos(x))
    testrtsafe(fdf, "sin(x) = 0"; x0=0.5, xl=-1.5, xh=1., root=0)

    fdf(x) = (sin(x), cos(x))
    testrtsafe(fdf, "sin(x) = 0"; x0=-1., xl=-1.5, xh=1., root=0)
end
