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
#include "PWGCF/FemtoUnited/DataModel/FemtoLambdasDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

#include "Framework/ASoA.h"

namespace o2::aod
{
namespace femtocascades
{
// columns for lambdas
DECLARE_SOA_COLUMN(XiMass, xiMass, float);       //! Mass of Lambda
DECLARE_SOA_COLUMN(OmegaMass, omegaMass, float); //! Mass of anti Lambda

// columns for Lambda bit masks
DECLARE_SOA_COLUMN(CascadeMask, cascadeMask, femtodatatypes::CascadeMaskType); //! Bitmask for Lambda selections

// columns for Lambda debug information
DECLARE_SOA_COLUMN(DauDCA, dauDCA, float);           //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Lambda transvers radius
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! x coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! y coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! z coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(KaonMass, kaonMass, float);       //! Lambda mass using Kaon hypothesis

// id columns for Lambda daughter tracks
DECLARE_SOA_INDEX_COLUMN_FULL(PosDauLambda, posDauLambda, int32_t, FUTracks, "_PosDauLambda"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegDauLambda, negDauLambda, int32_t, FUTracks, "_NegDauLambda"); //!

} // namespace femtocascades

// table for basic vzero information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FULambdas_001, "FUCASCADES", 1,
                                   o2::soa::Index<>,
                                   femtobase::CollisionId,
                                   femtobase::Pt,
                                   femtobase::Eta,
                                   femtobase::Phi,
                                   femtolambdas::LambdaMass,
                                   femtolambdas::AntiLambdaMass,
                                   femtolambdas::PosDauLambdaId,
                                   femtolambdas::NegDauLambdaId,
                                   femtobase::Theta<femtobase::Eta>,
                                   femtobase::Px<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Py<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                                   femtobase::P<femtobase::Pt, femtobase::Eta>);
using FULambdas = FULambdas_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FULambdaMasks_001, "FULAMBDAMASKS", 1,
                                   femtolambdas::LambdaMask);
using FULambdaMasks = FULambdaMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FULambdaExtras_001, "FULAMBDAEXTRAS", 1,
                                   femtolambdas::DauDCA,
                                   femtolambdas::DecayVtxX,
                                   femtolambdas::DecayVtxY,
                                   femtolambdas::DecayVtxZ,
                                   femtolambdas::TransRadius,
                                   femtolambdas::KaonMass);
using FULambdaExtras = FULambdaExtras_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_
