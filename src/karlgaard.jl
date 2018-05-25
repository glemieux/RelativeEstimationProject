using StaticArrays

# Karlgaard's formulation of the relative equations of motion
function sphericalkarlgaard(t, x::SVector, u::SVector, v::SVector, μ::Float64)

# Pull out parameters
ρ = x(1)
θ = x(3)
ϕ = x(5)
r = x(7)
f = x(9) # what is f?
ω = x(10)

dρ = x(2)
dθ = x(4)
dϕ = x(6)
dr = x(8)

# Define mu as gravparameter type?
Ga = [1 0 0; 0 1/ρ 0; 0 0 1/ρ]

# Common denominator
Ψ = (r^2 + ρ^2 + 2*r*ρ*cos(θ)*cos(ϕ))^(3/2)

# Spherical equations of motion
fρ = (ω^2 + 2*ω*dθ + dθ^2)*ρ*cos(ϕ)^2 + ρ*dϕ^2 - μ*(r*cos(θ)*cos(ϕ) + ρ)/Ψ + μ*cos(θ)*cos(ϕ)/r^2
fθ = 2*(ω + dθ)*dϕ*tan(ϕ) - dω - 2*(ω + dθ)*dρ/ρ + μ*(r*sin(θ)*sec(ϕ))/ρ*Ψ - μ*sin(θ)*sec(ϕ)/ρ*r^2
fϕ = -0.5*sin(2*ϕ)*(ω + dθ)^2 - 2*dρ*dϕ/ρ - μ*cos(θ)*sin(ϕ)/ρ*r^2 + μ*r*cos(θ)*sin(ϕ)/ρ*Ψ
fa = [fρ;fθ;fϕ]

# Need control noise term here

# Process noise term - Why did Jimmy multiply by gaussian random variables?  
v = diagm(randn(3))*sqrt.(v)

# Full relative motion acceleration with controls, η = {ρ;θ;ϕ}
ddη = fa + Ga*(u*v)
ddρ = ddη(1)
ddθ = ddη(2)
ddϕ = ddη(3)
ddr = r*ω^2 - μ/r^2
dω = -2*dr*ω/r

# State deriviative (velocity and acceleration)
dx = [dρ;ddρ;dθ;ddθ;dϕ;ddϕ;dr;ddr;ω;dω]
return dx

end
