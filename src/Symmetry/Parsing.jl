"Functions for parsing Crystals / SymOps from text / .cif files"

"Strips trailing uncertainty values from a String, then parses as a Float"
function _parse_cif_float(str::String) :: Float64
    i = findfirst('(', str)
    str = isnothing(i) ? str : str[1:i-1]
    return parse(Float64, str)
end


function parse_number_or_fraction(s)
    # Parse a number or fraction
    cnt = count(==('/'), s)
    if cnt == 0
        return parse(Float64, s)
    elseif cnt == 1
        n, d = parse.(Float64, split(s, '/'))
        return n/d
    else
        error("Cannot parse '$s' as a number.")
    end
end


function _parse_op(str::AbstractString) :: SymOp
    D = 3
    R = zeros(D, D)
    T = zeros(D)

    @assert length(split(str, ",")) == D

    # The substring s describes the i'th row of the output
    for (i, s) in enumerate(split(str, ","))
        # Remove all spaces
        s = replace(s, " " => "")
        # In weird circumstances, double negatives may be generated. Get rid of them.
        s = replace(s, "--" => "+")
        # Trick to split on either + or - symbols
        s = replace(s, "-" => "+-")

        # Each term t describes the numerical value for one element of R or T
        for t in split(s, '+', keepempty=false)
            j = findfirst(last(t), "xyz")
            if isnothing(j)
                T[i] += parse_number_or_fraction(t)
            else
                coeff = t[begin:end-1]
                R[i, j] += begin
                    if coeff == ""
                        1
                    elseif coeff == "-"
                        -1
                    else
                        parse_number_or_fraction(coeff)
                    end
                end
            end
        end
    end

    # Wrap translations
    T = mod.(T, 1)

    return SymOp(R, T)
end

"""
    Crystal(filename::AbstractString; symprec=1e-5)

Parse a `Crystal` from a `.cif` file located at the path of `filename`.
"""
function Crystal(filename::AbstractString; symprec=nothing)
    cif = Cif(Path(filename))
    # For now, assumes there is only one data collection per .cif
    cif = cif[first(keys(cif))]

    a = _parse_cif_float(cif["_cell_length_a"][1])
    b = _parse_cif_float(cif["_cell_length_b"][1])
    c = _parse_cif_float(cif["_cell_length_c"][1])
    α = _parse_cif_float(cif["_cell_angle_alpha"][1])
    β = _parse_cif_float(cif["_cell_angle_beta"][1])
    γ = _parse_cif_float(cif["_cell_angle_gamma"][1])
    lat_vecs = lattice_vectors(a, b, c, α, β, γ)

    geo_table = get_loop(cif, "_atom_site_fract_x")
    xs = map(_parse_cif_float, geo_table[!, "_atom_site_fract_x"])
    ys = map(_parse_cif_float, geo_table[!, "_atom_site_fract_y"])
    zs = map(_parse_cif_float, geo_table[!, "_atom_site_fract_z"])
    sitetypes = geo_table[!, "_atom_site_type_symbol"]
    sitetypes = Vector{String}(sitetypes)
    sitelabels = geo_table[!, "_atom_site_label"]
    unique_atoms = Vec3.(zip(xs, ys, zs))

    # Try to infer symprec from coordinate strings
    # TODO: Use uncertainty information if available from .cif
    if isnothing(symprec)
        strs = vcat(geo_table[!, "_atom_site_fract_x"], geo_table[!, "_atom_site_fract_y"], geo_table[!, "_atom_site_fract_z"])
        elems = vcat(xs, ys, zs)
        # guess fractional errors by assuming each elem is a fraction with simple denominator (2, 3, or 4)
        errs = map(elems) do x
            c = 12
            return abs(rem(x*c, 1, RoundNearest)) / c
        end
        (err, i) = findmax(errs)
        if err < 1e-12
            println("Precision parameter is unspecified, but all coordinates seem to be simple fractions.")
            println("Setting symprec=1e-12.")
            symprec = 1e-12
        elseif 1e-12 < err < 1e-4
            @printf "Precision parameter is unspecified, but coordinate string '%s' seems to have error %.1e.\n" strs[i] err
            symprec = 15err
            @printf "Setting symprec=%.1e.\n" symprec
        else
            println("Error: Please specify an explicit `symprec` parameter to load this file, '$filename'")
            return Nothing
        end
    end

    # By default, assume P1 symmetry
    symmetries = [ SymOp(I(3), @SVector zeros(3)) ]
    # There's two possible headers symmetry operations could be listed under
    for sym_header in ("_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz")
        if sym_header in keys(cif)
            sym_table = get_loop(cif, sym_header)
            symmetries = map(_parse_op, sym_table[!, sym_header])
        end
    end

    # By default, assume P1 symmetry
    groupnum = 1
    # Two possible headers for symmetry group number
    for group_header in ("_space_group_it_number", "_symmetry_int_tables_number")
        if group_header in keys(cif)
            groupnum = parse(Int, cif[group_header][1])
        end
    end

    hall_symbol = nothing
    hall_header = "_space_group_name_Hall"
    if hall_header in keys(cif)
        hall_symbol = cif[hall_header]
    end

    # Symmetry preferences: Explicit List > Hall Symbol > Infer
    if length(symmetries) > 1
        # Use explicitly provided symmetries
        return Crystal(lat_vecs, unique_atoms, sitetypes, symmetries; symprec)
    elseif !isnothing(hall_symbol)
        # Read symmetries from database for Hall symbol
        return Crystal(lat_vecs, unique_atoms, sitetypes, hall_symbol; symprec)
    else
        # Infer the symmetries automatically
        return Crystal(lat_vecs, unique_atoms, sitetypes; symprec)
    end
end