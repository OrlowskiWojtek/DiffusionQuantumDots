#include "include/visualiser.hpp"
// off clang format - order is important
// clang-format off
#include <morph/Visual.h>
#include <morph/VisualDataModel.h>
#include <morph/GridVisual.h>
#include <morph/Grid.h>
#include <morph/ColourBarVisual.h>
// clang-format on

void WalkersVisualiser::make_surf_plot(boost::multi_array<double, 2> &psi, int n) {
    morph::Visual v(800, 500, "Walkers distribution after simulation");

    v.showTitle(true);
    // Create a grid to show in the scene
    unsigned int Nside = n;
    constexpr morph::vec<float, 2> grid_spacing = {0.01f, 0.01f};

    morph::Grid grid(Nside, Nside, grid_spacing);

    std::vector<float> data(psi.num_elements(), 0.);
    std::transform(psi.data(), psi.data() + psi.num_elements(), data.begin(), [](double val){
        return static_cast<float>(val);
    });
    morph::vec<float, 3> offset = {-0.5f, -0.5f, 0.0f};

    morph::scale<float> zscale;
    zscale.do_autoscale = true;

    auto gv = std::make_unique<morph::GridVisual<float>>(&grid, offset);

    v.bindmodel(gv);
    gv->setSizeScale(1.);
    gv->gridVisMode = morph::GridVisMode::RectInterp;
    gv->zScale = zscale;
    gv->showborder (true);
    gv->border_thickness = 0.25f;
    gv->border_colour = morph::colour::black;
    gv->setScalarData(&data);
    gv->cm.setType(morph::ColourMapType::Twilight);
    gv->finalize();
    auto gvp = v.addVisualModel(gv);

    // Adding colorbar
    offset = {0.6f, -0.3f, 0.0f};
    auto cbv =  std::make_unique<morph::ColourBarVisual<float>>(offset);

    v.bindmodel (cbv);
    cbv->orientation = morph::colourbar_orientation::vertical;
    cbv->tickside = morph::colourbar_tickside::right_or_below;
    // Copy colourmap and scale to colourbar visual
    cbv->cm = gvp->cm;
    cbv->scale = gvp->colourScale;
    // Now build it
    cbv->finalize();
    v.addVisualModel (cbv);

    v.keepOpen();
}
