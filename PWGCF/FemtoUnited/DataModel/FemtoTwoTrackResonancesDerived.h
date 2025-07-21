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

/// \file FemtoTwoTrackResonancesDerived.h
/// \brief resonance tables
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTWOTRACKRESONANCESDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTWOTRACKRESONANCESDERIVED_H_

#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"

namespace o2::aod
{
namespace femtotwotrackresonances
{
// columns for resonance bit masks
DECLARE_SOA_COLUMN(MaskResonance, maskResonance, femtodatatypes::TwoTrackResonaceMaskType); //! Bitmask for resonance selections

DECLARE_SOA_COLUMN(MomentumPosDaughter, momentumPosDaughter, float); //! Bitmask for resonance selections
DECLARE_SOA_COLUMN(MomentumNegDaughter, momentumNegDaughter, float); //! Bitmask for resonance selections

DECLARE_SOA_INDEX_COLUMN_FULL(PosDauResonance, posDauResonance, int32_t, FUTracks, "_PosDauResonance"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegDauResonance, negDauResonance, int32_t, FUTracks, "_NegDauResonance"); //!

} // namespace femtotwotrackresonances

// table for phis
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUPhis_001, "FUPHIS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::CollisionId,
                                   femtobase::stored::Pt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtotwotrackresonances::PosDauResonanceId,
                                   femtotwotrackresonances::NegDauResonanceId,
                                   femtotwotrackresonances::MomentumPosDaughter,
                                   femtotwotrackresonances::MomentumNegDaughter,
                                   femtobase::dynamic::P<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Py<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Pz<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FUPhis = FUPhis_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUPhiMasks_001, "FUPHIMASKS", 1,
                                   femtotwotrackresonances::MaskResonance);
using FUPhiMasks = FUPhiMasks_001;

// table for kstars
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUKstars_001, "FUKSTARS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::CollisionId,
                                   femtobase::stored::SignedPt, // +1 for kstar and -1 for kstarbar
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtotwotrackresonances::PosDauResonanceId,
                                   femtotwotrackresonances::NegDauResonanceId,
                                   femtotwotrackresonances::MomentumPosDaughter,
                                   femtotwotrackresonances::MomentumNegDaughter,
                                   femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
                                   femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FUKstars = FUKstars_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUKstarMasks_001, "FUKSTARMASKS", 1,
                                   femtotwotrackresonances::MaskResonance);
using FUKstarMasks = FUKstarMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FURhos_001, "FURHOS", 1,
                                   o2::soa::Index<>,
                                   femtobase::stored::CollisionId,
                                   femtobase::stored::Pt,
                                   femtobase::stored::Eta,
                                   femtobase::stored::Phi,
                                   femtobase::stored::Mass,
                                   femtotwotrackresonances::PosDauResonanceId,
                                   femtotwotrackresonances::NegDauResonanceId,
                                   femtotwotrackresonances::MomentumPosDaughter,
                                   femtotwotrackresonances::MomentumNegDaughter,
                                   femtobase::dynamic::P<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Px<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Py<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Pz<femtobase::stored::Pt, femtobase::stored::Eta>,
                                   femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FURhos = FURhos_001;
DECLARE_SOA_TABLE_STAGED_VERSIONED(FURhoMasks_001, "FURHOMASKS", 1,
                                   femtotwotrackresonances::MaskResonance);
using FURhoMasks = FURhoMasks_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTWOTRACKRESONANCESDERIVED_H_
