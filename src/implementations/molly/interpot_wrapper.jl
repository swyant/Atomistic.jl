# -----------------------------------------------------------------------------
# Integration with InteratomicPotentials
# -----------------------------------------------------------------------------
const E_DIM = dimension(u"eV")
const EN_DIM = dimension(u"kJ * mol^-1") # need to account for annoying per atom case
const L_DIM = dimension(u"Å")

const DEF_EUNIT = u"eV"
const DEF_LUNIT = u"Å"

const NOUNIT_T = typeof(NoUnits)

struct InteratomicPotentialInter{P<:AbstractPotential}
    potential::P
    energy_units::Unitful.Unitlike 
    length_units::Unitful.Unitlike
    force_units::Unitful.Unitlike

    # internal constructor, ensuring energy units and length units have correct dimensions
    function InteratomicPotentialInter(
                                pot::P, 
                                eu::Union{NOUNIT_T, Unitful.Units{UE,E_DIM,nothing}, Unitful.Units{UE,EN_DIM,nothing}} = DEF_EUNIT,
                                lu::Union{NOUNIT_T, Unitful.Units{UL,L_DIM,nothing}} = DEF_LUNIT) where {P,UE,UL}
        if eu != NoUnits
            fu = eu/lu
        else
            fu = NoUnits
        end

        return new{P}(pot,eu,lu,fu) 
    end
end

function AtomsCalculators.forces( 
                      sys::AbstractSystem,
                      inter::InteratomicPotentialInter;
                      neighbors = nothing,
                      n_threads = Threads.nthreads())
    forces = InteratomicPotentials.force(sys,inter.potential)

    # initial profiling didn't show huge performance hit from unit conversion 
    # but! that may be because other parts of IP.jl are very slow, e.g. neighbor list construction
    if eltype(forces[1]) <: Unitful.Quantity 
        if inter.energy_units != NoUnits
            forces = [uconvert.(inter.force_units, fi)
                     for fi in forces]
        else
            forces = [ustrip.(fi)
                     for fi in forces]
        end
    elseif eltype(forces[1]) <: Real && inter.energy_units != NoUnits
        forces = forces * inter.force_units
    end

    forces
end

function AtomsCalculators.potential_energy( 
                                sys::AbstractSystem,
                                inter::InteratomicPotentialInter;
                                neighbors = nothing,
                                n_threads = Threads.nthreads())

    energy = InteratomicPotentials.potential_energy(sys,inter.potential)

    if typeof(energy) <: Unitful.Quantity 
        if inter.energy_units != NoUnits
            energy = uconvert(inter.energy_units, energy)
        else
            energy = ustrip(energy)
        end
    elseif typeof(energy) <: Real && inter.energy_units != NoUnits
        energy = energy * inter.energy_units
    end

    energy
end

AtomsCalculators.energy_unit(inter::InteratomicPotentialInter) = inter.energy_units
AtomsCalculators.length_unit(inter::InteratomicPotentialInter) = inter.length_units
