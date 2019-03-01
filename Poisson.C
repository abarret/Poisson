#include "Poisson.h"

#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"

using namespace tbox;
using namespace hier;
using namespace pdat;
using namespace geom;
namespace AMR
{
Poisson::Poisson(const std::string& object_name,
                 const tbox::Dimension& dim,
                 std::shared_ptr<tbox::Database> input_db,
                 std::shared_ptr<geom::CartesianGridGeometry> grid_geom)
{
    return;
}
void
Poisson::initializeLevelData(const std::shared_ptr<hier::PatchHierarchy>& /*hierarchy*/,
                             const int /*level_number*/,
                             const double /*init_data_time*/,
                             const bool /*can_be_refined*/,
                             const bool /*initial_time*/,
                             const std::shared_ptr<hier::PatchLevel>& /*old_level = nullptr*/,
                             const bool /*allocate_data = true*/)
{
    return;
}

void
Poisson::resetHierarchyConfiguration(const std::shared_ptr<hier::PatchHierarchy>& /*hierarchy*/,
                                     const int /*coarsest_level*/,
                                     const int /*finest_level*/)
{
    return;
}
}
