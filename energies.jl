#calculate energies
function kinetic_energy(m1, m2, w1, w2, pos1, pos2)

    # velocity of mass 1
    v1x = -w1 * pos1.y
    v1y =  w1 * pos1.x
    v1_squared = v1x^2 + v1y^2

    # relative position of mass 2 to mass 1
    rx = pos2.x - pos1.x
    ry = pos2.y - pos1.y

    # velocity of mass 2
    v2x = v1x - w2 * ry
    v2y = v1y + w2 * rx
    v2_squared = v2x^2 + v2y^2

    return 0.5*m1*v1_squared + 0.5*m2*v2_squared
end

function potential_energy(m1,m2,pos1,pos2,g) # mgh
    potential_energy = m1 * g * pos1.y + m2 * g * pos2.y
    return potential_energy
end

function total_energy(m1,m2,w1,w2,pos1,pos2,g)
    total_energy = kinetic_energy(m1,m2,w1,w2, pos1, pos2) + potential_energy(m1,m2,pos1,pos2,g)
    return total_energy
end

#apply energy calculations on list
function compute_energy_trajectory(positions1, positions2, w1s, w2s, m1, m2, g)

    kinetic_energies = []
    potential_energies = []
    total_energies = []
    
    for i in eachindex(positions1)
        ke = kinetic_energy(m1,m2,w1s[i],w2s[i],positions1[i],positions2[i])
        pe = potential_energy(m1,m2,positions1[i],positions2[i],g)
        te = ke + pe
        push!(kinetic_energies, ke)
        push!(potential_energies, pe)
        push!(total_energies, te)
    end
    
    return kinetic_energies, potential_energies, total_energies
end

#calculate standard deviations 
function standard_deviation(x)
    n = length(x)
    mu = sum(x) / n

    var = 0.0
    for xi in x
        var += (xi - mu)^2
    end

    var /= (n - 1)  

    return sqrt(var)
end
