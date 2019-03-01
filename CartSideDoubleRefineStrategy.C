// Filename: CartSideDoubleRefineStrategy.cpp
// Created on 30 Apr 2008 by Boyce Griffith
//
// Copyright (c) 2002-2017, Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

/////////////////////////////// INCLUDES /////////////////////////////////////

#include <ostream>
#include <set>
#include <vector>

#include "SAMRAI/hier/BoundaryBox.h"
#include "SAMRAI/hier/Box.h"
#include "SAMRAI/geom/CartesianPatchGeometry.h"
#include "SAMRAI/geom/CartesianSideDoubleConservativeLinearRefine.h"
#include "SAMRAI/hier/CoarseFineBoundary.h"
#include "SAMRAI/hier/ComponentSelector.h"
#include "SAMRAI/hier/Index.h"
#include "SAMRAI/hier/IntVector.h"
#include "SAMRAI/hier/Patch.h"
#include "SAMRAI/hier/PatchHierarchy.h"
#include "SAMRAI/hier/PatchLevel.h"
#include "SAMRAI/xfer/RefineAlgorithm.h"
#include "SAMRAI/hier/RefineOperator.h"
#include "SAMRAI/xfer/RefineSchedule.h"
#include "SAMRAI/pdat/SideData.h"
#include "SAMRAI/pdat/SideVariable.h"
#include "SAMRAI/hier/Variable.h"
#include "SAMRAI/hier/VariableContext.h"
#include "SAMRAI/hier/VariableDatabase.h"
#include "SAMRAI/tbox/Array.h"
#include "SAMRAI/tbox/Utilities.h"
#include "CartSideDoubleRefineStrategy.h"

using namespace SAMRAI;
using namespace tbox;
using namespace geom;
using namespace hier;
using namespace pdat;
using namespace xfer;

/////////////////////////////// NAMESPACE ////////////////////////////////////

/////////////////////////////// STATIC ///////////////////////////////////////

namespace
{
static const int REFINE_OP_STENCIL_WIDTH = 1;
static const int GHOST_WIDTH_TO_FILL = 1;
}

/////////////////////////////// PUBLIC ///////////////////////////////////////

CartSideDoubleRefineStrategy::CartSideDoubleRefineStrategy()
{
    return;
} // CartSideDoubleRefineStrategy

CartSideDoubleRefineStrategy::~CartSideDoubleRefineStrategy()
{
    return;
} // ~CartSideDoubleRefineStrategy

void
CartSideDoubleRefineStrategy::setPhysicalBoundaryConditions(Patch& /*patch*/,
                                                                      const double /*fill_time*/,
                                                                      const IntVector& /*ghost_width_to_fill*/)
{
    // intentionally blank
    return;
} // setPhysicalBoundaryConditions

IntVector
CartSideDoubleRefineStrategy::getRefineOpStencilWidth(const Dimension& dim) const
{
    return IntVector(dim, REFINE_OP_STENCIL_WIDTH);
} // getRefineOpStencilWidth

void
CartSideDoubleRefineStrategy::preprocessRefine(Patch& /*fine*/,
                                                         const Patch& /*coarse*/,
                                                         const Box& /*fine_box*/,
                                                         const IntVector& /*ratio*/)
{
    // intentionally blank
    return;
} // preprocessRefine

void
CartSideDoubleRefineStrategy::postprocessRefine(Patch& fine,
                                                          const Patch& coarse,
                                                          const Box& fine_box,
                                                          const IntVector& ratio)
{
    // intentionally blank
    return;
} // postprocessRefine

/////////////////////////////// PROTECTED ////////////////////////////////////

/////////////////////////////// PRIVATE //////////////////////////////////////
