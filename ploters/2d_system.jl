using GLMakie, CairoMakie
using DelimitedFiles


link_cameras_axis3(f; step=.01) = begin
  axes = filter(x -> x isa Axis3, f.content)

  for i ∈ 1:length(axes)
    lift(axes[i].azimuth, axes[i].elevation) do az, el
      for j ∈ 1:length(axes)
        i == j && continue
        if abs(az - axes[j].azimuth[]) > step
          axes[j].azimuth = az
        end
        if abs(el - axes[j].elevation[]) > step
          axes[j].elevation = el
        end
      end
    end
  end
  f
end
##

#filename = "../build/results/2d_harmonic_oscillator/1st_state_2d_anisotropic"
filename = "../build/evolution.dqmc.dat"
data = readdlm(filename, comments = true)

file = readlines(filename)
xmin = 0.
xmax = 0.
nbins = 0

for line in file
    if(contains(line, "xmin"))
        xmin = parse(Float64, split(line, ":")[2])
    end
    if(contains(line, "xmax"))
        xmax = parse(Float64, split(line, ":")[2])
    end
    if(contains(line, "nbins"))
        nbins = parse(Int64, split(line, ":")[2])
    end
end

x = y = LinRange(xmin, xmax, nbins)
# TODO: read from file
V(x,y) = 0.5 * (0.64 * (x^2) + 0.16 * (y^2))
Vmat = [V(xv, yv) for xv in x, yv in y]

##

GLMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis3(fig[1,1], 
               xlabel = "x [a.u]",
               ylabel = "y [a.u]")

    ax_v = Axis3(fig[1,1])
    hidedecorations!(ax_v)
    hidezdecorations!(ax)

    cm = surface!(ax, x, y, data,
             colormap = :coolwarm,
             label = "Ψ(x,y)")

    cmv = surface!(ax_v, x, y, Vmat,
             colormap = (:thermal, 0.5),
             transparency = true,
             label = "Potencjał harmoniczny")

    ticks_psi = [ round(val, digits = 3) for val in collect(LinRange(findmin(data)[1], findmax(data)[1], 6)) ]
    ticks_v =   [ round(val, digits = 3) for val in collect(LinRange(findmin(Vmat)[1], findmax(Vmat)[1], 6)) ]

    Colorbar(fig[2,1], cm, label = "Ψ(x,y)", vertical = false, ticks = ticks_psi[2:end-1])
    Colorbar(fig[3,1], cmv, label = "V(x,y)", vertical = false, ticks = ticks_v[2:end-1])

    xlims!(ax, (xmin, xmax))
    ylims!(ax, (xmin, xmax))
    zlims!(ax, (0, nothing))

    xlims!(ax_v, (xmin, xmax))
    ylims!(ax_v, (xmin, xmax))

    link_cameras_axis3(fig)
    
    #save("plots/2d_ground_state.pdf", fig)
    display(fig)
end
