module Hydraulics
include("code/PrismaticSection.jl")
include("code/CircularSection.jl")

# c = UniformPrismaticSection(5.0, 0.016, 2.0, 1.5, 1.5, 0.004)
# b = UniformSymmetricTrapezoidalSection(5.0, 0.016, 2.0, 1.5, 0.004)
# println(typeof(b))
# c = solve_uniform_flow(b)
# println(b.A)

c = UniformCircularSection(0.7, 0.014, 2.0, 0.0004)
println(c.h) # 0.68
end # module
