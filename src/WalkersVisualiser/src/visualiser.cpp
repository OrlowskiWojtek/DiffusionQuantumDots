#include "include/visualiser.hpp"
// off clang format - order is important
// clang-format off
#include <morph/Visual.h>
#include <morph/VisualDataModel.h>
#include <morph/GridVisual.h>
#include <morph/Grid.h>
#include <morph/ColourBarVisual.h>
// clang-format on

void WalkersVisualiser::make_surf_plot(boost::multi_array<double, 2> &psi,
                                       boost::multi_array<double, 2> &total_psi,
                                       int n) {
    morph::Visual v(800, 500, "Walkers distribution after simulation");

    v.showTitle(true);
    // Create a grid to show in the scene
    unsigned int Nside = n;
    constexpr morph::vec<float, 2> grid_spacing = {0.01f, 0.01f};

    morph::Grid grid(Nside, Nside, grid_spacing);

    std::vector<float> data(psi.num_elements(), 0.);
    std::transform(psi.data(), psi.data() + psi.num_elements(), data.begin(), [](double val) {
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
    gv->showborder(true);
    gv->border_thickness = 0.25f;
    gv->border_colour = morph::colour::black;
    gv->setScalarData(&data);
    gv->cm.setType(morph::ColourMapType::Twilight);
    gv->finalize();
    auto gvp = v.addVisualModel(gv);

    // Adding colorbar
    offset = {0.6f, -0.3f, 0.0f};
    auto cbv = std::make_unique<morph::ColourBarVisual<float>>(offset);

    v.bindmodel(cbv);
    cbv->orientation = morph::colourbar_orientation::vertical;
    cbv->tickside = morph::colourbar_tickside::right_or_below;
    // Copy colourmap and scale to colourbar visual
    cbv->cm = gvp->cm;
    cbv->scale = gvp->colourScale;
    // Now build it
    cbv->finalize();
    v.addVisualModel(cbv);

    // TODO: move to separate functions
    // total_grid
    morph::Grid total_grid(Nside, Nside, grid_spacing);

    std::vector<float> tot_data(total_psi.num_elements(), 0.);
    std::transform(total_psi.data(), total_psi.data() + total_psi.num_elements(), tot_data.begin(), [](double val) {
        return static_cast<float>(val);
    });
    morph::vec<float, 3> second_offset = {-0.5f, -1.5f, -0.f};

    morph::scale<float> second_zscale;
    second_zscale.do_autoscale = true;

    auto second_gv = std::make_unique<morph::GridVisual<float>>(&grid, second_offset);

    v.bindmodel(second_gv);
    second_gv->setSizeScale(1.);
    second_gv->gridVisMode = morph::GridVisMode::RectInterp;
    second_gv->zScale = zscale;
    second_gv->showborder(true);
    second_gv->border_thickness = 0.25f;
    second_gv->border_colour = morph::colour::black;
    second_gv->setScalarData(&tot_data);
    second_gv->cm.setType(morph::ColourMapType::Twilight);
    second_gv->finalize();
    auto sec_gvp = v.addVisualModel(second_gv);

    // Adding colorbar
    offset = {0.6f, -1.3f, 0.f};
    auto second_cbv = std::make_unique<morph::ColourBarVisual<float>>(offset);

    v.bindmodel(second_cbv);
    second_cbv->orientation = morph::colourbar_orientation::vertical;
    second_cbv->tickside = morph::colourbar_tickside::right_or_below;
    second_cbv->cm = sec_gvp->cm;
    second_cbv->scale = sec_gvp->colourScale;
    second_cbv->finalize();
    v.addVisualModel(second_cbv);

    v.keepOpen();
}
