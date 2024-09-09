#=
Tests for the staticAtoms method, as well as an implicit check that the Molly.system
conversion works appropriately.
=#
println(pwd())
init_sys_hfo2 = load_system("./implementations/molly/nonorthog_TiAL_example.xyz")
mod_sys_hfo2 = staticAtoms(init_sys_hfo2);

@test eltype(position(mod_sys_hfo2)) <: SVector
@test eltype(velocity(mod_sys_hfo2)) <: SVector

# TODO: stop spot-checking, check all values
@test position(mod_sys_hfo2)[1][2] == position(init_sys_hfo2)[1][2]
@test position(mod_sys_hfo2)[2][1] == position(init_sys_hfo2)[2][1]
@test position(mod_sys_hfo2)[2][3] == position(init_sys_hfo2)[2][3]
@test velocity(mod_sys_hfo2)[1][1] == velocity(mod_sys_hfo2)[1][1]

@test atomic_symbol(mod_sys_hfo2) == atomic_symbol(init_sys_hfo2)
@test atomic_mass(mod_sys_hfo2) == atomic_mass(init_sys_hfo2)
@test atomic_number(mod_sys_hfo2) == atomic_number(mod_sys_hfo2)

molly_sys_static = Molly.System(mod_sys_hfo2; energy_units=u"eV", force_units=u"eV/Å");
@test typeof(molly_sys_static) <: Molly.System
@test position(molly_sys_static)[1][1] == position(init_sys_hfo2)[1][1]
@test position(molly_sys_static)[1][3] == position(init_sys_hfo2)[1][3]
@test position(molly_sys_static)[2][1] == position(init_sys_hfo2)[2][1]

# This should work after https://github.com/JuliaMolSim/Molly.jl/pull/168
molly_sys_default = Molly.System(init_sys_hfo2; energy_units=u"eV", force_units=u"eV/Å")
@test typeof(molly_sys_default) <: Molly.System
@test position(molly_sys_default)[1][1] == position(init_sys_hfo2)[1][1]
@test position(molly_sys_default)[1][3] == position(init_sys_hfo2)[1][3]
@test position(molly_sys_default)[2][1] == position(init_sys_hfo2)[2][1]
