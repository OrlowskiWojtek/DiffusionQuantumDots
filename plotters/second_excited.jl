using CairoMakie
using DelimitedFiles

filepath = "../build/results/2nd_excited_nodes"
data = readdlm(filepath)

##

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "Położenie węzła", ylabel = "Końcowa wartość estymatora wzrostu")

    lines!(ax, data[:, 1], data[:, 2], color = :blue)
    scatter!(ax, data[:, 1], data[:, 2], color = :black)

    save("plots/2nd_excited_nodes.pdf", fig)
    display(fig)
end

##
filepath = "../build/results/2nd_excited_walkers"

x = Vector{Vector{Float64}}
x = LinRange(-5, 5, 200)
evolutions = Vector{Vector{Float64}}(undef, 0)
counter = 1
loading = false
temp_evo = Vector{Float64}(undef, 200)

for line in eachline(filepath)
    if(contains(line, "energy:"))
        loading = true
        temp_evo = Vector{Float64}(undef, 200)
        counter = 1
        continue
    end
    if(loading)
        temp_evo[counter] = parse(Float64, split(line, "\t")[2])
        counter += 1
        if(counter > 200)
            push!(evolutions, temp_evo)
            loading = false
        end
    end
end

##

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis3(fig[1,1], xlabel = "Pozycja węzłów", ylabel = "x", zlabel = "Końcowy rozkład wędrowców")

    nodes = data[:, 1]
    xs = x
    evo = [evolutions[i][j] for i in eachindex(nodes), j in eachindex(x)]
    
    cmap = cgrad(:bluesreds, length(nodes))
    for (idx, node) in enumerate(nodes)
        lines!(ax, [node], xs, evolutions[idx], color = cmap[idx], linewidth = 2)
    end

    lines!(ax, [1.07895], xs, evolutions[9], color = :black, linewidth = 2)

    save("plots/2nd_excited_walkers.pdf", fig)
    display(fig)
end
