using CairoMakie
using DelimitedFiles

#filepath = "../build/results/ground_energies"
filepath = "../build/averaged_energies.dqmc.dat"

data = readdlm(filepath, comments = true);
times = data[:, 1]
population = data[:, 2]
avg_mixed = data[:, 3]
avg_growth = data[:, 4]
mixed = data[:, 5]
growth = data[:, 6]
#

CairoMakie.activate!()

ntarget = split(readlines(filepath)[1], ":")[2]
#ntarget = 5000
#
with_theme(theme_latexfonts()) do
    fig = Figure(size = (1024, 768));
    col_mixed = :dodgerblue
    col_growth = :seagreen3
    col_walkers = :firebrick

    ax_growth = Axis(fig[1,1])
    ax_mixed = Axis(fig[3,1])
    ax_population = Axis(fig[5,1], xlabel = "Czas [a. u.]")

    lines!(ax_growth, times, growth, color = (col_growth, 0.5), label = "Estymator wzrostu [meV]")
    lines!(ax_growth, times, avg_growth, color = col_growth, label = "Średnia z estymatora wzrostu [meV]")

    lines!(ax_mixed, times, mixed, color = (col_mixed, 0.5), label = "Estymator mieszany [meV]")
    lines!(ax_mixed, times, avg_mixed, color = col_mixed, label = "Średnia z estymatora mieszanego [meV]")

    lines!(ax_population, times, population, color = (col_walkers, 0.5), label = "Liczba żywych wędrowców")
    hlines!(ax_population, [ntarget], color = col_walkers, label = "Celowana liczba wędrowców",  linestyle = :dash)

    Legend(fig[2,1], ax_growth, orientation = :horizontal, framevisible = false)
    Legend(fig[4,1], ax_mixed, orientation = :horizontal, framevisible = false)
    Legend(fig[6,1], ax_population, orientation = :horizontal, framevisible = false)

    linkxaxes!(ax_mixed, ax_growth) 
    linkxaxes!(ax_growth, ax_population) 

    #save("plots/1d_1el_ground_energies.pdf", fig)
    display(fig)
end
