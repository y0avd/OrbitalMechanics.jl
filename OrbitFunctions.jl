#=
Author: Yoav Dekel
=#

using LinearAlgebra

# Orbit Struct Functionality
struct Orbit
    a::Float64
    e::Float64
    ex::Float64
    ey::Float64
    i::Float64
    Ω::Float64
    θ::Float64
    ω::Float64
end

function eMagOrbit(a,e,i,Ω,θ,ω)
    ex = e*cos(ω)
    ey = e*sin(ω)

    return Orbit(a,e,ex,ey,i,Ω,θ,ω)
end

function eVecOrbit(a,ex,ey,i,Ω,θ)
    e = norm([ex,ey])
    ω = atan(ey,ex)

    return Orbit(a,e,ex,ey,i,Ω,θ,ω)
end

# Orbit Functions
function Orbit2Cart(μ,orbit::Orbit)
    ν = orbit.θ - orbit.ω
    p = orbit.a*(1-orbit.e^2)

    DCM = EA2DCM([orbit.Ω,orbit.i,orbit.ω],[3,1,3])

    r = [cos(ν),sin(ν),0]
    r = r*p/(1 + orbit.e*cos(ν))
    r = DCM'*r

    v = [-sin(ν),orbit.e + cos(ν),0]
    v = v*sqrt(μ/p)
    v = DCM'*v

    return [r;v]
end

function Cart2Orbit(μ,cart)
    r = cart[1:3]
    v = cart[4:6]

    h = cross(r,v)

    e = cross(v,h)/μ - r/norm(r)

    if e != 0
        ê = e/norm(e)
    else
        ê = [1 0 0]
    end

    a = (norm(h)^2/μ)/(1-norm(e)^2)

    ẑ = [0,0,1]

    ĥ = h/norm(h)

    println(ĥ)

    n = cross(ẑ,ĥ)

    if norm(n) != 0
        n̂ = n/norm(n)
    else
        n̂ = [1,0,0]
    end

    i = acos(h[3]/norm(h))
    Ω = atan(ĥ[1],-ĥ[2])

    cosω = ê[1]*n̂[1] + ê[2]*n̂[2]
    sinω = dot(cross(n̂,ê),ĥ)
    ω = atan(sinω,cosω)

    cosν = dot(ê,r)/norm(r)
    sinν = dot(cross(ê,r),ĥ)/norm(r)
    ν = atan(sinν,cosν)

    return eMagOrbit(a,norm(e),i,Ω,θ,ω)
end

function Kepler(M,e,ϵ = 1e-12,nᵢ = 1e3)
    # ϵ- tolerance
    # nᵢ- maximum iterations

    M = mod(M,2*pi)
    E = M

    for _ in 1:nᵢ
        f = M - E + e*sin(E)
        if abs(f) < ϵ
            cosν = (cos(E) - e)/(1 - e*cos(E))
            sinν = (sqrt(1 - e^2)*sin(E))/(1 - e*cos(E))
            return mod(atan(sinν,cosν),2pi)
        end

        Δf = -1 + e*cos(E)
        E = E - f/Δf
    end

    println("Error in Kepler(): Max iterations ($nᵢ) reached")
    return nothing # error
end

# Orbit Utils

function EA2DCM(EA,sequence)
    DCM = Matrix(1.0I,3,3)

    for i in 3:-1:1
        DCM = DCM*RFunction(EA[i],sequence[i])
    end

    return DCM
end

function RFunction(angle,axis)
    if axis == 1
        return [1 0 0;0 cos(angle) sin(angle);0 -sin(angle) cos(angle)]
    elseif axis == 2
        return [cos(angle) 0 -sin(angle);0 1 0;sin(angle) 0 cos(angle)]
    elseif axis == 3
        return [cos(angle) sin(angle) 0;-sin(angle) cos(angle) 0;0 0 1]
    end
end

# Propagators

"""
    This function is an Orbit Propagator for a spacecraft under J2 perturbations where the force model
        is using Cartesian coordinates

    This must be called with the ODEProblem() function, as it is using the ! syntax with du inside the function

    Inputs: (units are not important as long as consistent)
        du = vector that is derivative of u vector (relavent for calling in a DifferentialEquations.jl environment)
        u = vector with components: [r₁,r₂,r₃,v₁,v₂,v₃], the position and velocity vectors (in x y z Cartesian cordinates)
        consts = constants in the problem: [μ,J₂,R] which are the Gravitational parameter, J₂ pertubation constant and the radius of the planet
        t = time (relavent for calling in a DifferentialEquations.jl environment)
"""
function OrbitJ2PropCartesian!(du,u,consts::Tuple{Number,Number,Number},t)
    x,y,z = u
    μ,J₂,R = consts

    du[1:3] = u[4:6]

    r = sqrt(x^2 + y^2 + z^2)

    cc = (-3/2)μ*J₂*R^2/r^5 # placeholder

    du[4] = (-μ/r^3 + cc*(1 - 5z^2/r^2))x
    du[5] = (-μ/r^3 + cc*(1 - 5z^2/r^2))y
    du[6] = (-μ/r^3 + cc*(3 - 5z^2/r^2))z
end

"""
    This function is an Orbit Propagator for a spacecraft under J2 perturbations where the force model
        is using Cartesian coordinates

    This must be called with the ODEProblem() function, as it is using the ! syntax with du inside the function

    Inputs: (units are not important as long as consistent)
        du = vector that is derivative of u vector (relavent for calling in a DifferentialEquations.jl environment)
        u = vector with components: [a, e, i, Ω, ω, M], the position and velocity vectors (in x y z Cartesian cordinates)
        consts = constants in the problem: [μ,J₂,R] which are the Gravitational parameter, J₂ pertubation constant and the radius of the planet
        t = time (relavent for calling in a DifferentialEquations.jl environment)
"""
function OrbitJ2PropKeplerian!(du,u,consts,t)
    a,e,i,Ω,ω,M = u
    μ,J₂,R = consts
    
    ν = Kepler(M,e)

    DCM = EA2DCM([Ω,i,ν+ω],[3,1,3])

    p = a*(1 - e^2)
    b = a*sqrt(1 - e^2)
    h = sqrt(p*μ)
    r = p/(1 + e*cos(ν))

    f0 = [zeros(5,1);sqrt(μ/a^3)]

    B = (1/h)*[2e*(a^2)sin(ν) (2a^2)p/r 0;
        p*sin(ν) (p+r)cos(ν)+r*e 0;
        0 0 cos(ν+ω)r;
        0 0 sin(ν+ω)r/sin(i);
        -cos(ν)p/e (p+r)sin(ν)/e -sin(ν+ω)r/tan(i);
        b*cos(ν)p/(a*e)-2b*r/a -(p+r)sin(ν)b/(a*e) 0]

    x,y,z = DCM'*[r,0,0]

    aJ₂ = zeros(3,1) # malloc

    cc = (-3/2)μ*J₂*R^2/r^5 # placeholder

    aJ₂[1] = cc*(1 - 5z^2/r^2)x
    aJ₂[2] = cc*(1 - 5z^2/r^2)y
    aJ₂[3] = cc*(3 - 5z^2/r^2)z

    du[1:6] = f0 + B*(DCM*aJ₂)
end