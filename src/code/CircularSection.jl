export UniformCircularSection, CriticalCircularSection

include("Flow.jl")

mutable struct UniformCircularSection <: HydraulicSection
    @add_flow_fields
    @add_circular_geometry_fields
    So::Float64
    Sf::Float64
    function UniformCircularSection(Q::Float64, n::Float64, D::Float64, So::Float64)
        this = new()
        this.Q = Q
        this.n = n
        this.D = D
        this.So = this.Sf = So
        this.h = solve_uniform_flow(this)
        circular_section_hydraulic_properties(this, this.h)
        this
    end
end

mutable struct CriticalCircularSection <: HydraulicSection
    @add_flow_fields
    @add_circular_geometry_fields
    Sc::Float64
    function CriticalCircularSection(Q::Float64, n::Float64, D::Float64)
        this = new()
        this.Q = Q
        this.n = n
        this.D = D
        this.h = solve_critical_flow(this)
        this.Sc = ((Q * n) / (this.A * this.R^(2.0 / 3.0)))^2.0
        circular_section_hydraulic_properties(this, this.h)
        this
    end
end
