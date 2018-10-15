using Test, NewtonRaphsonSafe, Logging

global_logger(ConsoleLogger(stderr, Logging.Debug))

include("rtsafe.jl")
include("rtnewton.jl")
