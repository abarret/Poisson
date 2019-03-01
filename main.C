#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include <math.h>
#include <memory>
#include <stdio.h>
#include <stdlib.h>

#include "SAMRAI/SAMRAI_config.h"

// Headers for basic SAMRAI objects
#include "SAMRAI/appu/VisItDataWriter.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/math/HierarchyCellDataOpsReal.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/Dimension.h"
#include "SAMRAI/tbox/IOStream.h"
#include "SAMRAI/tbox/InputDatabase.h"
#include "SAMRAI/tbox/InputManager.h"
#include "SAMRAI/tbox/PIO.h"
#include "SAMRAI/tbox/SAMRAIManager.h"
#include "SAMRAI/tbox/SAMRAI_MPI.h"
#include "SAMRAI/tbox/Utilities.h"

// Headers for major algorithm/data structure objects
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/hier/CoarsenOperator.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/mesh/BergerRigoutsos.h"
#include "SAMRAI/mesh/ChopAndPackLoadBalancer.h"
#include "SAMRAI/mesh/GriddingAlgorithm.h"
#include "SAMRAI/mesh/StandardTagAndInitialize.h"
#include "SAMRAI/xfer/CoarsenAlgorithm.h"
#include "SAMRAI/xfer/CoarsenSchedule.h"

#include "CartSideDoubleRefineStrategy.h"
#include "Poisson.h"

using namespace SAMRAI;
using namespace tbox;
using namespace mesh;
using namespace geom;
using namespace appu;
using namespace hier;
using namespace pdat;
using namespace xfer;
using namespace AMR;
using std::string;

inline double
f(const std::vector<double>& x)
{
    return std::sin(2 * M_PI * x[0]) * std::sin(2 * M_PI * x[1]);
}
void fillU(std::shared_ptr<PatchHierarchy> patch_hierarchy, const int u_idx, const int ei);
void applyLaplacian(std::shared_ptr<PatchHierarchy> patch_hierarchy,
                    std::shared_ptr<SideVariable<double> > u_var,
                    const int u_idx,
                    const int r_idx);
int getCellWeightPatchDescriptorIndex(std::shared_ptr<PatchHierarchy> patch_hierarchy);
void printMatrix(std::shared_ptr<PatchHierarchy> patch_hierarchy, const int r_idx);

static std::ofstream mat_file;

int
main(int argc, char** argv)
{
    const Dimension dim(2);
    SAMRAI_MPI::init(&argc, &argv);
    SAMRAIManager::initialize();
    SAMRAIManager::startup();
    string input_filename;
    if (argc != 2)
    {
        tbox::pout << "USAGE: " << argv[0] << " <input filename> "
                   << "[options]\n\n";
        tbox::SAMRAI_MPI::abort();
    }
    else
    {
        input_filename = argv[1];
    }
    plog << "input_filename = " << input_filename << "\n";
    std::shared_ptr<MemoryDatabase> input_db = InputManager::getManager()->parseInputFile(input_filename);
    auto main_db = std::static_pointer_cast<MemoryDatabase>(input_db->getDatabase("Main"));

    string log_file_name = "Poisson.log";
    if (main_db->keyExists("log_file_name"))
    {
        log_file_name = main_db->getString("log_file_name");
    }

    string visit_dump_dirname;
    int visit_number_procs_per_file = 1;
    if (main_db->keyExists("visit_dump_dirname"))
    {
        visit_dump_dirname = main_db->getString("visit_dump_dirname");
    }
    if (main_db->keyExists("visit_number_procs_per_file"))
    {
        visit_number_procs_per_file = main_db->getInteger("visit_number_procs_per_file");
    }

    auto grid_geometry =
        std::make_shared<CartesianGridGeometry>(dim, "CartesianGeometry", input_db->getDatabase("CartesianGeometry"));

    auto patch_hierarchy =
        std::make_shared<PatchHierarchy>("PatchHierarchy", grid_geometry, input_db->getDatabase("PatchHierarchy"));

    auto visit_data_writer = std::make_shared<VisItDataWriter>(
        dim, "Oldroyd B Visit Writer", visit_dump_dirname, visit_number_procs_per_file);

    auto poisson = std::make_shared<Poisson>("Poisson", dim, nullptr, grid_geometry);

    auto tag_and_init_ops = std::make_shared<StandardTagAndInitialize>(
        "StandardTagAndInitialize", poisson.get(), input_db->getDatabase("StandardTagAndInitialize"));

    // Gridding algorithm stuff...
    auto box_generator = std::make_shared<BergerRigoutsos>(dim);

    auto load_balancer = std::make_shared<ChopAndPackLoadBalancer>(
        dim, "ChopAndPackLoadBalancer", input_db->getDatabase("LoadBalancer"));

    auto gridding_algorithm = std::make_shared<GriddingAlgorithm>(patch_hierarchy,
                                                                  "GriddingAlgorithm",
                                                                  input_db->getDatabase("GriddingAlgorithm"),
                                                                  tag_and_init_ops,
                                                                  box_generator,
                                                                  load_balancer);

    // Print database and variable database contents to log filename
    plog << "\nCheck input data and variables before simulation:\n";
    plog << "Input database...\n";
    input_db->printClassData(plog);
    plog << "\nVariable Database...\n";
    VariableDatabase::getDatabase()->printClassData(plog);

    // Initialize Data on patches
    std::vector<int> tag_buffer_array(patch_hierarchy->getMaxNumberOfLevels());
    for (int i = 0; i < patch_hierarchy->getMaxNumberOfLevels(); ++i)
    {
        tag_buffer_array[i] = 1;
        pout << "i = " << i << " tag_buffer = " << tag_buffer_array[i] << "\n";
    }

    double loop_time = 0.0;

    gridding_algorithm->makeCoarsestLevel(loop_time);

    bool done = false;
    bool initial_time = true;
    for (int ln = 0; patch_hierarchy->levelCanBeRefined(ln) && !done; ln++)
    {
        gridding_algorithm->makeFinerLevel(tag_buffer_array[ln], initial_time, 0, loop_time);
        done = !(patch_hierarchy->finerLevelExists(ln));
    }
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    auto directions = IntVector::getZero(dim);
    directions(0) = 1;
    auto u_var = std::make_shared<SideVariable<double> >(dim, "U", directions);
    auto u_draw_var = std::make_shared<CellVariable<double> >(dim, "U Draw");
    auto r_var = std::make_shared<SideVariable<double> >(dim, "Grad U", directions);

    std::shared_ptr<VariableContext> ctx = var_db->getContext("U::context");

    const int u_idx = var_db->registerVariableAndContext(u_var, ctx, IntVector::getOne(dim));
    const int u_draw_idx = var_db->registerVariableAndContext(u_draw_var, ctx, IntVector::getZero(dim));
    const int r_idx = var_db->registerVariableAndContext(r_var, ctx, IntVector::getZero(dim));

    visit_data_writer->registerPlotQuantity(u_draw_var->getName(), "SCALAR", u_draw_idx);

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        std::shared_ptr<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        level->allocatePatchData(u_idx);
        level->allocatePatchData(u_draw_idx);
        level->allocatePatchData(r_idx);
    }

    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        std::shared_ptr<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            std::shared_ptr<Patch> patch = *p;
            const Box& box = patch->getBox();
            auto u_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(u_idx));
            auto pgeom = std::static_pointer_cast<CartesianPatchGeometry>(patch->getPatchGeometry());
            const double* dx = pgeom->getDx();
            const double* xlow = pgeom->getXLower();
            const IntVector& idx_low = box.lower();

            SideIterator cend(SideGeometry::end(box, 0));
            for (SideIterator c(SideGeometry::begin(box, 0)); c != cend; ++c)
            {
                SideIndex idx = *c;
                std::vector<double> x(dim.getValue());
                for (int d = 0; d < dim.getValue(); ++d)
                    x[d] = xlow[d] + dx[d] * (idx(d) - idx_low(d) + d == 0 ? 0 : 0.5);
                (*u_data)(idx) = f(x);
            }
        }
    }

    // Loop through and count degrees of freedom
    int dof = 0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        std::shared_ptr<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            std::shared_ptr<Patch> patch = *p;
            const Box& box = patch->getBox();
            SideIterator send(SideGeometry::end(box, 0));
            for (SideIterator s(SideGeometry::begin(box, 0)); s != send; ++s)
            {
                dof++;
            }
        }
    }
    pout << dof << "\n";
    mat_file.open("matrix_vals");
    // Loop through e_i and apply laplacian
    for (int i = 0; i < dof; ++i)
    {
        pout << "Creating and filling row: " << i << "\n";
        fillU(patch_hierarchy, u_idx, i);
        applyLaplacian(patch_hierarchy, u_var, u_idx, r_idx);
        printMatrix(patch_hierarchy, r_idx);
    }
    // Print out result -> this is column of matrix

    visit_data_writer->writePlotData(patch_hierarchy, 0, 0.0);
    tbox::pout << "\n\n##########################################\n";
    tbox::pout << " Printing out visualization files...\n";
    tbox::pout << "\n###########################################\n";

    SAMRAIManager::shutdown();
    SAMRAIManager::finalize();
    SAMRAI_MPI::finalize();

    return 0;
}

void
printMatrix(std::shared_ptr<PatchHierarchy> patch_hierarchy, const int r_idx)
{
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        std::shared_ptr<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            std::shared_ptr<Patch> patch = *p;
            const Box& box = patch->getBox();
            auto r_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(r_idx));
            SideIterator send(SideGeometry::end(box, 0));
            for (SideIterator s(SideGeometry::begin(box, 0)); s != send; ++s)
            {
                const SideIndex& idx = *s;
                double val = (*r_data)(idx);
                mat_file << std::to_string(val) << " ";
            }
        }
    }
    mat_file << "\n";
    return;
}

void
fillU(std::shared_ptr<PatchHierarchy> patch_hierarchy, const int u_idx, const int ei)
{
    int i = 0;
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        std::shared_ptr<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            std::shared_ptr<Patch> patch = *p;
            const Box& box = patch->getBox();
            auto u_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(u_idx));
            u_data->fillAll(0.0);
            SideIterator send(SideGeometry::end(box, 0));
            for (SideIterator si(SideGeometry::begin(box, 0)); si != send; ++si)
            {
                const SideIndex& idx = *si;
                if (i == ei) (*u_data)(idx) = 1.0;
                i++;
            }
        }
    }
    return;
}

void
applyLaplacian(std::shared_ptr<PatchHierarchy> patch_hierarchy,
               std::shared_ptr<SideVariable<double> > u_var,
               const int u_idx,
               const int r_idx)
{
    std::shared_ptr<CartesianGridGeometry> grid_geometry =
        std::static_pointer_cast<CartesianGridGeometry>(patch_hierarchy->getGridGeometry());
    const Dimension& dim = patch_hierarchy->getDim();
    // Fill in ghost cells
    if (patch_hierarchy->getFinestLevelNumber() > 0)
    {
        int fill_level_num = 1;
        int coarser_level = 0;
        std::shared_ptr<CoarsenOperator> coarsen_op(
            grid_geometry->lookupCoarsenOperator(u_var, "CONSERVATIVE_COARSEN"));
        CoarsenAlgorithm coarsen_alg(dim);
        coarsen_alg.registerCoarsen(u_idx, u_idx, coarsen_op, IntVector::getZero(dim));
        std::shared_ptr<CoarsenSchedule> coarsen_sched;
        coarsen_sched = coarsen_alg.createSchedule(patch_hierarchy->getPatchLevel(coarser_level),
                                                   patch_hierarchy->getPatchLevel(fill_level_num));
        coarsen_sched->coarsenData();
        std::shared_ptr<RefineOperator> refine_op(
            grid_geometry->lookupRefineOperator(u_var, "CONSERVATIVE_LINEAR_REFINE"));
        RefineAlgorithm ghost_fill_alg;
        ghost_fill_alg.registerRefine(u_idx, u_idx, u_idx, refine_op);
        std::shared_ptr<RefineSchedule> fill_ghosts_sched;
        auto patch_strategy = std::make_shared<CartSideDoubleRefineStrategy>();
        fill_ghosts_sched = ghost_fill_alg.createSchedule(
            patch_hierarchy->getPatchLevel(fill_level_num), coarser_level, patch_hierarchy, nullptr /*patch_strategy*/);
        fill_ghosts_sched->fillData(0.0, false);
    }
    // Do differencing
    for (int ln = 0; ln <= patch_hierarchy->getFinestLevelNumber(); ++ln)
    {
        std::shared_ptr<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        for (auto p = level->begin(); p != level->end(); ++p)
        {
            std::shared_ptr<Patch> patch = *p;
            const Box& box = patch->getBox();
            auto u_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(u_idx));
            auto r_data = std::static_pointer_cast<SideData<double> >(patch->getPatchData(r_idx));
            auto pgeom = std::static_pointer_cast<CartesianPatchGeometry>(patch->getPatchGeometry());
            const double* const dx = pgeom->getDx();
            SideIterator send(SideGeometry::end(box, 0));
            for (SideIterator si(SideGeometry::begin(box, 0)); si != send; ++si)
            {
                const SideIndex& idx = *si;
                IntVector t = IntVector::getZero(dim), r = IntVector::getZero(dim);
                t(1) = 1;
                r(0) = 1;
                (*r_data)(idx) = ((*u_data)(idx + t) + (*u_data)(idx - t) + (*u_data)(idx + r) + (*u_data)(idx - r) -
                                  4 * (*u_data)(idx)) /
                                 (dx[0] * dx[1]);
            }
        }
    }
    return;
}

int
getCellWeightPatchDescriptorIndex(std::shared_ptr<PatchHierarchy> patch_hierarchy)
{
    const Dimension& dim = patch_hierarchy->getDim();
    VariableDatabase* var_db = VariableDatabase::getDatabase();
    auto cc_var = std::make_shared<CellVariable<double> >(dim, "cc_var", 1);
    std::shared_ptr<VariableContext> cxt = var_db->getContext("WEIGHTS");
    int wgt_cc_idx = var_db->registerVariableAndContext(cc_var, cxt, IntVector(dim, 0));
    for (int ln = patch_hierarchy->getFinestLevelNumber(); ln >= 0; ln--)
    {
        std::shared_ptr<PatchLevel> level = patch_hierarchy->getPatchLevel(ln);
        if (!level->checkAllocated(wgt_cc_idx)) level->allocatePatchData(wgt_cc_idx);
        for (PatchLevel::Iterator p = level->begin(); p != level->end(); p++)
        {
            std::shared_ptr<Patch> patch = *p;
            const Box& patch_box = patch->getBox();
            auto pgeom = std::static_pointer_cast<CartesianPatchGeometry>(patch->getPatchGeometry());

            const double* dx = pgeom->getDx();
            double cell_vol = 1.0;
            for (int d = 0; d < dim.getValue(); ++d) cell_vol *= dx[d];
            auto wgt_cc_data = std::static_pointer_cast<CellData<double> >(patch->getPatchData(wgt_cc_idx));
            wgt_cc_data->fillAll(cell_vol);
        }
    }
    return wgt_cc_idx;
}
