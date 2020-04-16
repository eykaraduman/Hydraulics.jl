macro def(name, definition)
  return quote
      macro $(esc(name))()
          esc($(Expr(:quote, definition)))
      end
  end
end

@def add_flow_fields begin
    Q::Float64
    n::Float64
    h::Float64
    A::Float64
    P::Float64
    T::Float64
    R::Float64
    V::Float64
    hv::Float64
    Fr::Float64
    E::Float64
end

@def add_geometry_fields begin
    b::Float64
    m1::Float64
    m2::Float64
end

@def add_circular_geometry_fields begin
    D::Float64
    beta::Float64
end

global const g = 9.80665

abstract type HydraulicSection end

function prismatic_section_hydraulic_properties(o::T, y::Float64) where {T<:HydraulicSection}
     o.A = y * (o.b + y * (o.m1 + o.m2) / 2.0)
     o.P = o.b + y * (sqrt(1.0 + o.m1^2) + sqrt(1.0 + o.m2^2))
     o.R = o.A / o.P
     o.T = o.b + y * (o.m1 + o.m2)
     o.V = o.Q / o.A
     o.hv = o.V^2 / (2.0 * g)
     o.E = y + o.hv
     o.Fr = o.V / sqrt(g * o.A / o.T)
     (o.Q, o.A, o.R, o.So, o.n, o.T)
end

function circular_section_hydraulic_properties(o::T, y::Float64) where {T<:HydraulicSection}
     o.beta = 2.0 * acos(1.0 - 2.0 * y / o.D)
     o.A = o.D^2.0 * (o.beta - sin(o.beta)) / 8.0
     o.T = o.D * sin(o.beta / 2.0)
     o.P = o.D * o.beta / 2.0
     o.R = o.A / o.P
     o.V = o.Q / o.A
     o.hv = o.V^2 / (2.0 * g)
     o.E = y + o.hv
     o.Fr = o.V / sqrt(g * o.A / o.T)
     (o.Q, o.A, o.R, o.So, o.n, o.T)
end

function uniform_flow_height(o::T, y::Float64) where {T<:HydraulicSection}
    Q, A, R, So, n = typeof(o) == UniformPrismaticSection ?
    prismatic_section_hydraulic_properties(o, y) : circular_section_hydraulic_properties(o, y)
    Q - A * R ^(2.0 / 3.0) * sqrt(So) / n
end

function critical_flow_height(o::T, y::Float64) where {T<:HydraulicSection}
       Q, A, R, So, n, Tt = typeof(o) == UniformPrismaticSection ?
       prismatic_section_hydraulic_properties(o, y) : circular_section_hydraulic_properties(o, y)
       sqrt(g * A^3.0 / Tt) - Q
end

function calculate_prismatic_section_discharge(o::T, y::Float64) where {T<:HydraulicSection}
    A = y * (o.b + y * (o.m1 + o.m2) / 2.0)
    P = o.b + y * (sqrt(1.0 + o.m1^2) + sqrt(1.0 + o.m2^2))
    R = A / P
    A * R ^(2.0 / 3.0) * sqrt(o.So) / o.n
end

using Roots
# using ForwardDiff

function solve_uniform_flow(o::T) where {T<:HydraulicSection}
    x1 = typeof(o) == UniformPrismaticSection ? 10.00 : 0.95 * o.D
    f = y -> uniform_flow_height(o, y)
    fzero(f, 0.01, x1)
end

function solve_critical_flow(o::T) where {T<:HydraulicSection}
    x1 = typeof(o) == UniformPrismaticSection ? 10.00 : 0.95 * o.D
    f = y -> critical_flow_height(o, y)
    fzero(f, 0.01, x1)
end
