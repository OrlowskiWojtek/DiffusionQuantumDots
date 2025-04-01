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

filename = "../build/results/2d_harmonic_oscillator/1st_state_2d_anisotropic"
data = readdlm(filename, comments = true)

lines = readlines(filename)
xmin = 0.
xmax = 0.
nbins = 0

for line in lines
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

V(x,y) = 0.5 * (0.64 * (x^2) + 0.16 * (y^2))
Vmat = [V(xv, yv) for xv in x, yv in y]

##

GLMakie.activate!()

with_theme(theme_latexfonts()) do
    fig = Figure();
    ax = Axis3(fig[1,1], 
               xlabel = "x [a.u]",
               ylabel = "y [a.u]",
               zlabel = "Ψ(x,y)")

    ax_v = Axis3(fig[1,1])
    hidedecorations!(ax_v)

    surface!(ax, x, y, data,
             colormap = :bluesreds)

    surface!(ax_v, x, y, Vmat,
             colormap = (:acton, 0.6),
             transparency = true)
    link_cameras_axis3(fig)
    
    #save("2d_ground_state.pdf", fig)
    display(fig)
end

