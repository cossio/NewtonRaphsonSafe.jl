function testrtnewton(fdf, eqn; x0, root)
    print("solving " * eqn * " from x0 = $x0 (real root = $root) ... ")
    call_count = 0
    F(x) = begin call_count += 1; fdf(x) end
    x = rtnewton(F, x0)
    print("took $call_count function calls\n")
    @test x ≈ root  atol=1e-15
end


@testset "rtnewton" begin
    fdf(x) = (exp(x) - 3, exp(x))
    testrtnewton(fdf, "exp(x) = 1"; x0=5, root=log(3))
    testrtnewton(fdf, "exp(x) = 1"; x0=1, root=log(3))
    testrtnewton(fdf, "exp(x) = 1"; x0=-1, root=log(3))
end