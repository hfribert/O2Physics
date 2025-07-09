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

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOV0SDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOV0SDERIVED_H_

#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"

namespace o2::aod
{
namespace femtov0s
{
// columns for lambdas and K0s
DECLARE_SOA_COLUMN(LambdaMass, lambdaMass, float);         //! Mass of Lambda
DECLARE_SOA_COLUMN(AntiLambdaMass, antiLambdaMass, float); //! Mass of anti Lambda
DECLARE_SOA_COLUMN(KaonMass, kaonMass, float);             //! Lambda mass using Kaon hypothesis

// columns for bit masks
DECLARE_SOA_COLUMN(LambdaMask, lambdaMask, femtodatatypes::V0MaskType); //! Bitmask for Lambda selections
DECLARE_SOA_COLUMN(KaonMask, kaonMask, femtodatatypes::V0MaskType);     //! Bitmask for Lambda selections

// columns for debug information
DECLARE_SOA_COLUMN(CosPa, cosPa, float);             //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(DauDCA, dauDca, float);           //! Lambda daughter DCA at decay vertex
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Lambda transvers radius
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! x coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! y coordinate of Lambda decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! z coordinate of Lambda decay vertex

// id columns for Lambda daughter tracks
DECLARE_SOA_INDEX_COLUMN_FULL(PosDau, posDau, int32_t, FUTracks, "_PosDau"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegDau, negDau, int32_t, FUTracks, "_NegDau"); //!

} // namespace femtov0s

// table for basic lambda information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FULambdas_001, "FULAMBDAS", 1,
                                   o2::soa::Index<>,
                                   femtobase::CollisionId,
                                   femtobase::Pt,
                                   femtobase::Eta,
                                   femtobase::Phi,
                                   femtov0s::LambdaMass,
                                   femtov0s::AntiLambdaMass,
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::Theta<femtobase::Eta>,
                                   femtobase::Px<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Py<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                                   femtobase::P<femtobase::Pt, femtobase::Eta>);
using FULambdas = FULambdas_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FULambdaMasks_001, "FULAMBDAMASKS", 1,
                                   femtov0s::LambdaMask);
using FULambdaMasks = FULambdaMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FULambdaExtras_001, "FULAMBDAEXTRAS", 1,
                                   femtov0s::CosPa,
                                   femtov0s::DauDCA,
                                   femtov0s::DecayVtxX,
                                   femtov0s::DecayVtxY,
                                   femtov0s::DecayVtxZ,
                                   femtov0s::TransRadius,
                                   femtov0s::KaonMass);

using FULambdaExtras = FULambdaExtras_001;

// table for basic lambda information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUKshorts_001, "FUKSHORTS", 1,
                                   o2::soa::Index<>,
                                   femtobase::CollisionId,
                                   femtobase::Pt,
                                   femtobase::Eta,
                                   femtobase::Phi,
                                   femtov0s::KaonMass,
                                   femtov0s::PosDauId,
                                   femtov0s::NegDauId,
                                   femtobase::Theta<femtobase::Eta>,
                                   femtobase::Px<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Py<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                                   femtobase::P<femtobase::Pt, femtobase::Eta>);
using FUKshorts = FUKshorts_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUKshortMasks_001, "FUKSHORTMASKS", 1,
                                   femtov0s::KaonMask);
using FUKshortMasks = FUKshortMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUKshortExtras_001, "FUSHORTEXTRAS", 1,
                                   femtov0s::CosPa,
                                   femtov0s::DauDCA,
                                   femtov0s::DecayVtxX,
                                   femtov0s::DecayVtxY,
                                   femtov0s::DecayVtxZ,
                                   femtov0s::TransRadius);

using FUKshortExtras = FUKshortExtras_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOV0SDERIVED_H_
