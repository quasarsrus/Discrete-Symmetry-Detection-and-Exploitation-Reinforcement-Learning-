using Pkg
using Plots
using ViscousStreaming
using HDF5
using JLD2

 Re = 40
 ϵ = 0.1
 Ω = 1.0 # frequency (keep this equal to 1)
 Tp = 2π/Ω # one period of oscillation
 p = StreamingParams(ϵ,Re)
 #s = StreamingAnalytical(p)

 τ = 0.1 # Stokes number, should be small
 β = 0.95 # Density parameter. Less than 1 means heavier than fluid.
 p_inert = InertialParameters(tau=τ,beta=β,epsilon=ϵ,Re=Re)
 Ω = 1.0
 Tp = 2π/Ω
 Tmax = 50*Tp

 Δx = 0.02
 xlim = (-2.8,2.8)
 ylim = (-2.8,2.8)
 n = 75
 body = Circle(0.2,n)

 bl = BodyList()
 bL1 = deepcopy(body)
 bL2 = deepcopy(body)
 bR1 = deepcopy(body)
 bR2 = deepcopy(body)

"""
#BigSquare
# left 1 cylinder
cent = (-5.0,5.0)
 α = 0.0
TL1 = RigidTransform(cent,α)
TL1(bL1) # transform the body to the current configuration

# left 2 cylinder
cent = (-5.0,-5.0)
 α = 0.0
TL2 = RigidTransform(cent,α)
TL2(bL2) # transform the body to the current configuration

# right 1 cylinder
cent = (5.0,5.0)
 α = 0.0
TR1 = RigidTransform(cent,α)
TR1(bR1) # transform the body to the current configuration

# right 2 cylinder
cent = (5.0,-5.0)
 α = 0.0
TR2 = RigidTransform(cent,α)
TR2(bR2) # transform the body to the current configuration

push!(bl,bL1);
push!(bl,bL2);
push!(bl,bR1);
push!(bl,bR2);
"""

#Square
# left 1 cylinder
cent = (-2.0,2.0)
 α = 0.0
TL1 = RigidTransform(cent,α)
TL1(bL1) # transform the body to the current configuration

# left 2 cylinder
cent = (-2.0,-2.0)
 α = 0.0
TL2 = RigidTransform(cent,α)
TL2(bL2) # transform the body to the current configuration

# right 1 cylinder
cent = (2.0,2.0)
 α = 0.0
TR1 = RigidTransform(cent,α)
TR1(bR1) # transform the body to the current configuration

# right 2 cylinder
cent = (2.0,-2.0)
 α = 0.0
TR2 = RigidTransform(cent,α)
TR2(bR2) # transform the body to the current configuration

push!(bl,bL1);
push!(bl,bL2);
push!(bl,bR1);
push!(bl,bR2);

"""
#Trapezium
# left 1 cylinder
cent = (-1.0,2.0)
 α = 0.0
TL1 = RigidTransform(cent,α)
TL1(bL1) # transform the body to the current configuration

# left 2 cylinder
cent = (-2.0,-2.0)
 α = 0.0
TL2 = RigidTransform(cent,α)
TL2(bL2) # transform the body to the current configuration

# right 1 cylinder
cent = (1.0,2.0)
 α = 0.0
TR1 = RigidTransform(cent,α)
TR1(bR1) # transform the body to the current configuration

# right 2 cylinder
cent = (2.0,-2.0)
 α = 0.0
TR2 = RigidTransform(cent,α)
TR2(bR2) # transform the body to the current configuration

push!(bl,bL1);
push!(bl,bL2);
push!(bl,bR1);
push!(bl,bR2);
"""

"""
#Quadrilateral
# left 1 cylinder
cent = (-2.0,1.3)
 α = 0.0
TL1 = RigidTransform(cent,α)
TL1(bL1) # transform the body to the current configuration

# left 2 cylinder
cent = (-2.3,-1.0)
 α = 0.0
TL2 = RigidTransform(cent,α)
TL2(bL2) # transform the body to the current configuration

# right 1 cylinder
cent = (1.5,2.2)
 α = 0.0
TR1 = RigidTransform(cent,α)
TR1(bR1) # transform the body to the current configuration

# right 2 cylinder
cent = (0.1,-2.5)
 α = 0.0
TR2 = RigidTransform(cent,α)
TR2(bR2) # transform the body to the current configuration

push!(bl,bL1);
push!(bl,bL2);
push!(bl,bR1);
push!(bl,bR2);
"""

@time solver2 = FrequencyStreaming(Re,ϵ,Δx,xlim,ylim,bl);

ampvec = [ComplexF64[0.0,0.0] for i in 1:length(bl)];
 ampvec[1] = [1,0];  a2 = deepcopy(ampvec);
 ampvec[2] = [1,0];  a3 = deepcopy(ampvec); 
 ampvec[3] = [1,0];  a4 = deepcopy(ampvec); 
 ampvec[4] = [1,0];  a5 = deepcopy(ampvec);
 ampvec[1] = [0,0];  a6 = deepcopy(ampvec);
 ampvec[2] = [0,0];  a7 = deepcopy(ampvec);
 ampvec[3] = [0,0];  a8 = deepcopy(ampvec);
 ampvec[2] = [1,0];  a9 = deepcopy(ampvec);
 ampvec[1] = [1,0];  a10 = deepcopy(ampvec);
 ampvec[2] = [0,0];  a11 = deepcopy(ampvec);
 ampvec[3] = [1,0];  a12 = deepcopy(ampvec);
 ampvec[4] = [0,0];  a13 = deepcopy(ampvec);
 ampvec[1] = [0,0];  a14 = deepcopy(ampvec);
 ampvec[2] = [1,0];  a15 = deepcopy(ampvec);
 ampvec[3] = [0,0];  a16 = deepcopy(ampvec);

 actionspace = [a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16];

strdVuxy = []
strdVvxy = []
for i = 1: length(actionspace)
    soln = solver2(actionspace[i],bl);
    isoln = inertial_velocity(soln,p_inert);
    v̄L = lagrangian_mean_velocity(isoln);
    v̄Luxy, v̄Lvxy = interpolatable_field(v̄L,isoln.g);
    push!(strdVuxy,v̄Luxy)
    push!(strdVvxy,v̄Lvxy)
end

save("strdVuxy_sq.jld2", "data", strdVuxy)
save("strdVvxy_s.jld2", "data", strdVvxy)
