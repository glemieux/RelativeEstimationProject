using StaticArrays

# Karlgaard's formulation of the relative equations of motion
function sphericalkarlgaard(t, x::SVector, u::SVector, v::SVector, mu::Float64)
# Define mu as gravparameter type?

Ga = [1 0 0; 0 1/x(1) 0; 0 0 1/x(1)]

# Pull out parameters
ρ = x(1)
θ = x(3)
ϕ = x(5)
r = x(7)
f = x(9)

dρ = x(2)
dθ = x(4)
dϕ = x(6)
dr = x(8)
df = x(10)
ddf = x(11)



end
