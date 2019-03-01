#include "Poisson.h"

#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/hier/VariableContext.h"

using namespace tbox;
using namespace hier;
using namespace pdat;
using namespace geom;
namespace AMR
{
Poisson::Poisson(const std::string& object_name,
                 const Dimension& dim,
                 std::shared_ptr<tbox::Database> input_db,
                 std::shared_ptr<CartesianGridGeometry> grid_geom)
:       d_dim(dim), d_grid_geom(grid_geom)
{
/*    VariableDatabase* var_db = VariableDatabase::getDatabase();
    std::shared_ptr<VariableContext> ctx = var_db->getContext("Context");
    auto directions = IntVector::getZero(dim);
    directions(0) = 1;
    d_u_var = std::make_shared<SideVariable<double>>(d_dim, "U", directions);
    d_r_var = std::make_shared<SideVariable<double>>(d_dim, "Grad U", directions);
    
    d_u_idx = var_db->registerVariableAndContext(d_u_var, ctx, IntVector::getOne(d_dim));
    d_r_idx = var_db->registerVariableAndContext(d_r_var, ctx, IntVector::getZero(d_dim));*/
    return;
}
Poisson::~Poisson()
{
    return;
}

void Poisson::initializeLevelData(const std::shared_ptr<hier::PatchHierarchy>& /*hierarchy*/,
                             const int /*level_number*/,
                             const double /*init_data_time*/,
                             const bool /*can_be_refined*/,
                             const bool /*initial_time*/,
                             const std::shared_ptr<hier::PatchLevel>& /*old_level = nullptr*/,
                             const bool /*allocate_data = true*/)
{
    return;
}

void Poisson::resetHierarchyConfiguration(const std::shared_ptr< hier::PatchHierarchy >& /*hierarchy*/,
                                          const int /*coarsest_level*/,
                                          const int /*finest_level*/)
{
    return;
}

void Poisson::computeLaplace(const std::shared_ptr<PatchHierarchy>& hierarchy)
{
    
}

}
