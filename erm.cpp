/*
 * Implements 'Subbarrel Patterns in Somatosensory Cortical Barrels Can Emerge from
 * Local Dynamic Instabilities', Ermentrout, Simons and Land, PLOS Comp Biol 2009
 */

#define FLT float

#include <cstdint>
#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <utility>
#include <ctime>
#include <filesystem>

import sm.vec;
import sm.vvec;
import sm.algo;
import sm.config;
import sm.interval;
import mplot.colourmap;
import mplot.visual;
import mplot.hexgridvisual;
import mplot.graphvisual;

import rd_erm;

int main (int argc, char **argv)
{
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " /path/to/params.json [/path/to/logdir]" << std::endl;
        return 1;
    }
    std::string paramsfile (argv[1]);

    // Get parameters from config file
    sm::config conf(paramsfile);
    if (!conf.ready) {
        std::cerr << "Failed to read config " << paramsfile << ". Exiting.\n";
        return 1;
    }

    const std::uint32_t steps = conf.get<std::uint32_t> ("steps", 1000);
    if (steps == 0) {
        std::cerr << "Finished simulating 0 steps. Exiting.\n";
        return 1;
    }
    const std::uint32_t logevery = conf.get<std::uint32_t> ("logevery", 100);
    if (logevery == 0) {
        std::cerr << "Can't log every 0 steps. Exiting.\n";
        return 1;
    }
    bool overwrite_logs = conf.get<bool> ("overwrite_logs", false);
    std::string logpath = conf.get<std::string> ("logpath", "logs/erm2");
    if (argc == 3) {
        std::string argpath(argv[2]);
        std::cerr << "Overriding the config-given logpath " << logpath << " with " << argpath << std::endl;
        logpath = argpath;
        if (overwrite_logs == true) {
            std::cerr << "WARNING: You set a command line log path.\n"
                      << "       : Note that the parameters config permits the program to OVERWRITE LOG\n"
                      << "       : FILES on each run (\"overwrite_logs\" is set to true).\n";
        }
    }

    std::cout << "steps to simulate: " << steps << std::endl;

    const std::uint32_t plotevery = conf.get<std::uint32_t> ("plotevery", 10);
    const std::uint32_t win_width = conf.get<std::uint32_t> ("win_width", 340);
    std::uint32_t win_height = static_cast<std::uint32_t>(0.8824f * (float)win_width);

    // Instantiate the model object
    sm::rd_erm<FLT> RD;
    RD.svgpath = conf.get<std::string> ("svgpath", ""); // We'll do an elliptical boundary if svgpath is empty
    RD.ellipse_a = conf.get<FLT> ("ellipse_a", FLT{0.8});
    RD.ellipse_b = conf.get<FLT> ("ellipse_b", FLT{0.6});
    RD.logpath = logpath;
    RD.set_dt(static_cast<FLT>(conf.get<FLT> ("dt", FLT{0.00001})));
    RD.hextohex_d = conf.get<FLT> ("hextohex_d", FLT{0.01});
    RD.boundaryFalloffDist = conf.get<FLT> ("boundaryFalloffDist", FLT{0.01});
    RD.N = 1;
    RD.Dn = conf.get<FLT> ("Dn", FLT{0.3});
    RD.Dc = conf.get<FLT> ("Dc", FLT{0.3} * RD.Dn); // 0.3 * Dn
    RD.beta = conf.get<FLT> ("beta", FLT{5});
    RD.a = conf.get<FLT> ("a", FLT{1});
    RD.b = conf.get<FLT> ("b", FLT{1});
    RD.mu = conf.get<FLT> ("mu", FLT{1});
    RD.chi = conf.get<FLT> ("chi", FLT{-1});
    if (RD.chi < 0.0) {
        // Set chi to Dn, as in the paper (see linear analysis)
        RD.chi = RD.Dn; // or: conf.get<FLT> ("chi", 0.3);
    }

    // Allocate and initialise the model
    RD.hexspan = 3.0f  * (RD.ellipse_a > RD.ellipse_b ? RD.ellipse_a : RD.ellipse_b);
    std::cout << "RD.hexspan: " << RD.hexspan << ", Dn=" << RD.Dn << ", Dc=" << RD.Dc << ", chi=" << RD.chi << std::endl;
    RD.allocate();
    RD.init();

    // Autoscale the colour map on each step?
    bool do_autoscale = conf.get<bool> ("autoscale", false);

    // Now create a log directory if necessary, and exit on any failures.
    if (std::filesystem::is_directory (logpath) == false) {
        std::filesystem::create_directories (logpath);
        if (std::filesystem::is_directory (logpath) == false) {
            std::cerr << "Failed to create the logpath directory "
                      << logpath << " which does not exist.\n";
            return 1;
        }
    } else {
        // Directory DOES exist. See if it contains a previous run and
        // exit without overwriting to avoid confusion.
        if (overwrite_logs == false) {
            if ((std::filesystem::exists (logpath + "/params.json")
                 && std::filesystem::is_regular_file (logpath + "/params.json"))
                ||
                (std::filesystem::exists (logpath + "/positions.h5")
                 && std::filesystem::is_regular_file (logpath + "/positions.h5"))) {
                std::cerr << "Seems like a previous simulation was logged in " << logpath << ".\n"
                          << "Please clean it out manually, choose another directory or set\n"
                          << "overwrite_logs to true in your parameters config JSON file.\n";
                return 1;
            }
        }
    }

    // As RD.allocate() as been called (and log directory has been
    // created/verified ready), positions can be saved to file.
    RD.savePositions();

    // Set up the mplot::Visual object
    mplot::Visual<> v (win_width, win_height, "Ermentrout (Keller-Segel)");
    // You can set the default scene x/y/z offsets
    v.setSceneTrans (conf.get<float> ("x_default", 0.0f),
                     conf.get<float> ("y_default", 0.0f),
                     conf.get<float> ("z_default", -5.0f));

    // A 2D scaling to apply to the visuals
    float myscale = conf.get<float> ("size_scale", 1.0f);

    // A z position to place the hexgrid visuals
    float _Z = 0.0f;

    sm::vec<float, 3> spatOff = {0,0,_Z}; // spatial offset
    // Data scaling parameters
    float _m = 0.2;
    float _c = 0.0;
    sm::scale<FLT, float> cscale;
    cscale.set_params (_m, _c);

    // Set up a 3D map of the surface RD.n[0] using a mplot::HexGridVisual
    spatOff[0] -= myscale * (0.6 * (RD.hg->width()));
    auto hgv1 = std::make_unique<mplot::HexGridVisual<FLT>> (RD.hg.get(), spatOff);
    hgv1->set_parent (v.get_id());
    hgv1->scaleViewMatrix (myscale);
    hgv1->setScalarData (&RD.n[0]);
    // You can directly set VisualDataModel::zScale and ::colourScale:
    hgv1->zScale.set_params (_m/10.0f, _c/10.0f);
    // ...or use setters to copy one in:
    hgv1->setCScale (cscale);
    hgv1->cm.setType (mplot::ColourMapType::Jet);
    hgv1->hexVisMode = mplot::HexVisMode::Triangles;
    hgv1->addLabel ("n (axon density)", {-0.6f, RD.hg->width() / 1.9f, 0}, mplot::TextFeatures (0.12f, 64));
    hgv1->finalize();
    auto hgv1p = v.addVisualModel (hgv1);

    // Set up a 3D map of the surface RD.c[0]
    spatOff[0] *= -1;
    auto hgv2 = std::make_unique<mplot::HexGridVisual<FLT>> (RD.hg.get(), spatOff);
    hgv2->set_parent (v.get_id());
    hgv2->scaleViewMatrix (myscale);
    hgv2->setScalarData (&RD.c[0]);
    hgv2->zScale.set_params (_m/10.0f, _c/10.0f);
    hgv2->setCScale (cscale);
    hgv2->cm.setType (mplot::ColourMapType::Jet);
    hgv2->hexVisMode = mplot::HexVisMode::HexInterp;
    hgv2->addLabel ("c (chemoattractant)", {-0.7f, RD.hg->width() / 1.9f, 0}, mplot::TextFeatures (0.12f, 64));
    hgv2->finalize();
    auto hgv2p = v.addVisualModel (hgv2);

    // Set up a 2D graph with mplot::GraphVisual
    spatOff = {0.5f, -2.0f, 0.0f};
    auto graph1 = std::make_unique<mplot::GraphVisual<FLT>> (spatOff);
    graph1->set_parent (v.get_id());
    graph1->setlimits (0, steps * RD.get_dt(), 0, conf.get<float>("graph_ymax", 40000.0f));
    graph1->policy = mplot::stylepolicy::lines;
    graph1->ylabel = "Sum";
    graph1->xlabel = "Sim time (s)";
    graph1->prepdata ("n");
    graph1->prepdata ("c");
    graph1->finalize();
    auto graph1p = v.addVisualModel (graph1);

    // Set up the render clock
    std::chrono::steady_clock::time_point lastrender = std::chrono::steady_clock::now();

    // Start the loop
    bool doing = true;
    static constexpr bool debug_ranges = false;

    // Use a pointer cast to get access to the range() function in vvec with the std::vectors n[0] and c[0]
    auto n0 = static_cast<sm::vvec<FLT>*>(&RD.n[0]);
    auto c0 = static_cast<sm::vvec<FLT>*>(&RD.c[0]);

    while (doing && !v.readyToFinish()) {
        // Step the model:
        RD.step();

        if ((RD.stepCount % plotevery) == 0) {
            // Plot n and c
            if (do_autoscale == true) {
                hgv1p->colourScale.compute_scaling (n0->range());
                hgv2p->colourScale.compute_scaling (c0->range());
            }

            hgv1p->updateData (n0);
            hgv2p->updateData (c0);
            // Append to the 2D graph of sums:
            graph1p->append ((float)RD.stepCount * RD.get_dt(), RD.sum_n[0], 0);
            graph1p->append ((float)RD.stepCount * RD.get_dt(), RD.sum_c[0], 1);
        }

         // Save data every 'logevery' steps
        if ((RD.stepCount % logevery) == 0) {
            std::cout << "Logging data at step " << RD.stepCount << std::endl;
            RD.saveState();
        }

        // After a while, stop:
        if (RD.stepCount > steps) { doing = false; }

        // render the graphics not on every loop, but at about 60 Hz. On non-rendering loops, we only step the RD system.
        std::chrono::steady_clock::duration sincerender = std::chrono::steady_clock::now() - lastrender;
        if (std::chrono::duration_cast<std::chrono::milliseconds>(sincerender).count() > 17) { // 17 is about 60 Hz
            v.poll();
            v.render();
            lastrender = std::chrono::steady_clock::now();
        }
    }

    // Before exit, save data
    RD.saveState();

    // Add simulation runtime information to the config, before saving it out as params.json
    conf.set ("float_width", static_cast<std::uint32_t>(sizeof(FLT)));
    std::time_t curtime = std::time (nullptr);
    std::string tnow = std::asctime(std::localtime (&curtime));
    conf.set ("sim_ran_at_time", tnow.substr(0,tnow.size()-1));
    conf.set ("final_step", RD.stepCount);
    conf.set ("hextohex_d", RD.hextohex_d);
    conf.set ("dt", RD.get_dt());
    // Store the binary name and command argument into root, too.
    if (argc > 0) { conf.set("argv0", argv[0]); }
    if (argc > 1) { conf.set("argv1", argv[1]); }
    const std::string paramsCopy = logpath + "/params.json";
    conf.write (paramsCopy);
    if (conf.ready == false) {
        std::cerr << "Warning: Something went wrong writing a copy of the params.json: "
                  << conf.emsg << std::endl;
    }

    std::cout << "Finished simulating. Press Ctrl-q in graphics window to exit.\n";
    v.keepOpen();

    return 0;
}
