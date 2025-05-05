using CairoMakie
using DelimitedFiles

function read_dist(filename::String)
    walker_dist = Vector{Vector{Float64}}(undef, 0)
    times = Float64[]
    filelines = readlines(filename)
    xmin = 0.
    xmax = 0.
    nbins = 0
    for line in filelines
        if(contains(line, "xmin"))
            xmin = parse(Float64, split(line, ":")[2])
        end
        if(contains(line, "xmax"))
            xmax = parse(Float64, split(line, ":")[2])
        end
        if(contains(line, "nbins"))
            nbins = parse(Int64, split(line, ":")[2])
        end
        if(!contains(line, "#"))
            break
        end
    end
    
    temp_dist = Float64[]
    for line in filelines
        if(contains(line, "#"))
            if(contains(line, "time:"))
                push!(times, parse(Float64, split(line, ":")[2]))
            end
            continue
        end
        if(line == "")
            push!(walker_dist, temp_dist)
            temp_dist = Float64[]
            continue
        end
        push!(temp_dist, parse(Float64, split(line, "\t")[1]))
    end

    x = collect(LinRange(xmin, xmax, nbins))

    return x, walker_dist, times
end

function psi_0(x)
    return @. (0.8 / π)^(0.25) * exp(-0.8 * x^2 / 2.)
end

filename = "../build/evolution.dqmc.dat"
x, walker_dist, times = read_dist(filename)
##
#
with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis(fig[1,1:2], xlabel = "położenie [a.u.]", ylabel = "|Ψ(x)|²",
              ylabelcolor = :blue, title = "Stan podstawowy oscylatora harmonicznego.");

    cmap = cgrad(:managua, length(times), categorical = true)
    for (idx, walkers) in enumerate(walker_dist)
        lines!(ax, x, walkers, label = "czas = $(times[idx]) a.u.", color = cmap[idx])
    end

    lines!(ax, x, psi_0(x).^2,
           label = "Rozwiązanie dokładne",
           color = :red,
           linestyle = :dash)

    Legend(fig[1,3], ax, framevisible = false)

    display(fig)
end
