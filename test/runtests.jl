using Test, NewtonRaphsonSafe, Logging

const _LOGGER = ConsoleLogger(stderr, Logging.Debug)
global_logger(_LOGGER)

include("rtsafe.jl")
include("rtnewton.jl")
