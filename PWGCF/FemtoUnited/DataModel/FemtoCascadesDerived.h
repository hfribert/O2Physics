// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.femtobaseder
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file FemtoLambdasDerived.h
/// \brief v0 tables tables
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_

#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoV0sDerived.h"

#include "Framework/ASoA.h"

namespace o2::aod
{
namespace femtocascades
{
// columns for xis and omegas
DECLARE_SOA_COLUMN(XiMass, xiMass, float);       //! Mass of Lambda
DECLARE_SOA_COLUMN(OmegaMass, omegaMass, float); //! Mass of anti Lambda

// columns for Lambda bit masks
DECLARE_SOA_COLUMN(CascadeMask, cascadeMask, femtodatatypes::CascadeMaskType); //! Bitmask for Lambda selections

// columns for cascad debug information
DECLARE_SOA_COLUMN(CascadeCosPa, cascadeCosPa, float);             //! cosine of the poiting angle at decay vertex
DECLARE_SOA_COLUMN(CascadeDauDca, cascadeDauDca, float);           //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(CascadeTransRadius, cascadeTransRadius, float); //! Lambda transvers radius
DECLARE_SOA_COLUMN(LambdaCosPa, lambdaCosPa, float);               //! cosine of the poiting angle at decay vertex
DECLARE_SOA_COLUMN(LambdaDauDca, lambdaDauDca, float);             //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(LambdaTransRadius, lambdaTransRadius, float);   //! Lambda transvers radius
DECLARE_SOA_COLUMN(LambdaDcaToPv, lambdaDcaToPv, float);           //! Lambda transvers radius

// id columns for daughter lambda
// DECLARE_SOA_INDEX_COLUMN_FULL(LambdaDaughter, lambdaDaughter, int32_t, FULambdas, "_LambdaDaughter"); //!

// id columns for bachelor
DECLARE_SOA_INDEX_COLUMN_FULL(CascadeBachelor, cascadeBachelor, int32_t, FUTracks, "_CascadeBachelor"); //!
// id columns for Lambda daughter tracks
DECLARE_SOA_INDEX_COLUMN_FULL(PosDauLambda, posDauLambda, int32_t, FUTracks, "_PosDauLambda"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegDauLambda, negDauLambda, int32_t, FUTracks, "_NegDauLambda"); //!

} // namespace femtocascades

// table for basic vzero information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUCascades_001, "FUCASCADES", 1,
                                   o2::soa::Index<>,
                                   femtobase::CollisionId,
                                   femtobase::SignedPt,
                                   femtobase::Eta,
                                   femtobase::Phi,
                                   femtocascades::XiMass,
                                   femtocascades::OmegaMass,
                                   femtocascades::CascadeBachelorId,
                                   femtocascades::PosDauLambdaId,
                                   femtocascades::NegDauLambdaId,
                                   femtotracks::Pt<femtobase::SignedPt>,
                                   femtotracks::Sign<femtobase::SignedPt>,
                                   femtobase::P<femtobase::SignedPt, femtobase::Eta>,
                                   femtobase::Px<femtobase::SignedPt, femtobase::Eta>,
                                   femtobase::Py<femtobase::SignedPt, femtobase::Eta>,
                                   femtobase::Pz<femtobase::SignedPt, femtobase::Eta>,
                                   femtobase::Theta<femtobase::Eta>, );
using FUCascades = FUCascades_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUCascadeMasks_001, "FUCASCADEMASKS", 1,
                                   femtocascades::CascadeMask);
using FUCascadeMasks = FUCascadeMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUCascadeExtras_001, "FUCASCADEEXTRAS", 1,
                                   femtocascades::CascadeCosPa,
                                   femtocascades::CascadeDauDca,
                                   femtocascades::CascadeTransRadius,
                                   femtocascades::LambdaCosPa,
                                   femtocascades::LambdaDauDca,
                                   femtocascades::LambdaTransRadius,
                                   femtocascades::LambdaDcaToPv);
using FUCascadeExtras = FUCascadeExtras_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_
