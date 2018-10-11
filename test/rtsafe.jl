using LambertWExp


function testwrap(fdf, eqn; x0, xl, xh, root)
    print("solving " * eqn * " from x0 = $x0 ... ")
    call_count = 0
    F(x) = begin call_count += 1; fdf(x) end
    @test rtsafe(F, x0, xl, xh) ≈ root
    print("took $call_count function calls\n")
end


@testset "rtsafe" begin
    fdf(x) = begin W = lambertwexp(x); (W - 1, W / (1 + W)) end
    testwrap(fdf, "lambertwexp(x) = 1"; x0=5, xl=-10, xh=10, root=1)

    fdf(x) = (sin(x), cos(x))
    testwrap(fdf, "sin(x) = 0"; x0=0.5, xl=-1, xh=1, root=0)
end