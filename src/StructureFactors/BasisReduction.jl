function phase_averaged_elements(data, k::Vec3, sf::StructureFactor{N, NCorr, NAtoms}, ffdata::Vector{FormFactor{FFType}}) where {FFType, N, NCorr, NAtoms}
    elems = zero(SVector{NCorr, ComplexF64})
    knorm = norm(k)
    ffs = ntuple(i -> compute_form(knorm, ffdata[i]), NAtoms)
    rs = ntuple(i -> sf.crystal.latvecs * sf.crystal.positions[i], NAtoms)

    for j in 1:NAtoms, i in 1:NAtoms
        phase = exp(im*(k ⋅ (rs[j] - rs[i])))
        elems += phase * ffs[i] * ffs[j] * view(data, :, i, j)  # This view allocates
    end

    return elems
end
