export UniformPrismaticSection, CriticalPrismaticSection, GvfPrismaticSection

include("Flow.jl")

mutable struct UniformPrismaticSection <: HydraulicSection
    @add_flow_fields
    @add_geometry_fields
    So::Float64
    Sf::Float64
    function UniformPrismaticSection(Q::Float64, n::Float64, b::Float64, m1::Float64, m2::Float64, So::Float64)
        this = new()
        this.Q = Q
        this.n = n
        this.b = b
        this.m1 = m1
        this.m2 = m2
        this.So = this.Sf = So
        this.h = solve_uniform_flow(this)
        prismatic_section_hydraulic_properties(this, this.h)
        this
    end

    function UniformPrismaticSection(Q::Float64, n::Float64, b::Float64, m::Float64, So::Float64)
        this = UniformPrismaticSection(Q, n, b, m, m, So)
        this
    end

    function UniformPrismaticSection(Q::Float64, n::Float64, b::Float64, So::Float64)
        this = UniformPrismaticSection(Q, n, b, 0.0, 0.0, So)
        this
    end
end

mutable struct CriticalPrismaticSection <: HydraulicSection
    @add_flow_fields
    @add_geometry_fields
    Sc::Float64
    function CriticalPrismaticSection(Q::Float64, n::Float64, b::Float64, m1::Float64, m2::Float64)
        this = new()
        this.Q = Q
        this.n = n
        this.b = b
        this.m1 = m1
        this.m2 = m2
        this.h = solve_critical_flow(this)
        prismatic_section_hydraulic_properties(this, this.h)
        this.Sc = ((Q * n) / (this.A * this.R^(2.0 / 3.0)))^2.0
        this
    end

    function CriticalPrismaticSection(Q::Float64, n::Float64, b::Float64, m::Float64)
        this = CriticalPrismaticSection(Q, n, b, m, m)
        this
    end

    function CriticalPrismaticSection(Q::Float64, n::Float64, b::Float64)
        this = CriticalPrismaticSection(Q, n, b, 0.0, 0.0)
        this
    end
end

mutable struct GvfPrismaticSection <: HydraulicSection
    @add_flow_fields
    @add_geometry_fields
    So::Float64
    Sf::Float64
    function GvfPrismaticSection(Q::Float64, h::Float64, n::Float64, b::Float64, m1::Float64, m2::Float64, So::Float64)
        this = new()
        this.Q = Q
        this.h = h
        this.n = n
        this.b = b
        this.m1 = m1
        this.m2 = m2
        this.So = So
        this.Sf = ((Q * n) / (this.A * this.R^(2.0 / 3.0)))^2.0
        prismatic_section_hydraulic_properties(this, h)
        this
    end

    function GvfPrismaticSection(Q::Float64, h::Float64, n::Float64, b::Float64, m::Float64, So::Float64)
        this = GvfPrismaticSection(Q, h, n, m, m, So)
        this
    end

    function GvfPrismaticSection(Q::Float64, h::Float64, n::Float64, b::Float64, So::Float64)
        this = GvfPrismaticSection(Q, h, n, 0.0, 0.0, So)
        this
    end
end
