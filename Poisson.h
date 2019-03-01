#ifndef included_Poisson
#define included_Poisson

#include "SAMRAI/SAMRAI_config.h"

#include "SAMRAI/appu/BoundaryUtilityStrategy.h"
#include "SAMRAI/geom/CartesianGridGeometry.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/mesh/StandardTagAndInitStrategy.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/tbox/Database.h"
#include "SAMRAI/tbox/MessageStream.h"
#include "SAMRAI/tbox/Serializable.h"

#include <memory>
#include <vector>

using namespace SAMRAI;
namespace AMR
{
/**
  * Dummy implementation of StandardTagAndInitStrategy. Provides empty implementations of intializeLevelData and
 * resetHierarchyConfiguration.
  */
class Poisson : public mesh::StandardTagAndInitStrategy
{
public:
    Poisson(const std::string& object_name,
            const tbox::Dimension& dim,
            std::shared_ptr<tbox::Database> input_db,
            std::shared_ptr<geom::CartesianGridGeometry> grid_geom);
    
    ~Poisson();
    
    void initializeLevelData(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                             const int level_number,
                             const double init_data_time,
                             const bool can_be_refined,
                             const bool initial_time,
                             const std::shared_ptr<hier::PatchLevel>& old_level = nullptr,
                             const bool allocate_data = true) override;
    void resetHierarchyConfiguration(const std::shared_ptr<hier::PatchHierarchy>& hierarchy,
                                     const int coarsest_level,
                                     const int finest_level) override;

private:
};
}
#endif
