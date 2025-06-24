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

// columns for resonances
DECLARE_SOA_COLUMN(ResonanceMass, resonanceMass, float);                                //! Mass of resonance
DECLARE_SOA_COLUMN(ResonanceType, resonanceType, femtodatatypes::TwoTrackResonaceType); //! resonance type

// columns for resonance bit masks
DECLARE_SOA_COLUMN(ResonanceMask, resonanceMask, femtodatatypes::TwoTrackResonaceMaskType); //! Bitmask for resonance selections

DECLARE_SOA_INDEX_COLUMN_FULL(PosDauResonance, posDauResonance, int32_t, FUTracks, "_PosDauResonance"); //!
DECLARE_SOA_INDEX_COLUMN_FULL(NegDauResonance, negDauResonance, int32_t, FUTracks, "_NegDauResonance"); //!

} // namespace femtotwotrackresonances

// table for resonances
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUResos_001, "FURESOS", 1,
                                   o2::soa::Index<>,
                                   femtobase::CollisionId,
                                   femtobase::Pt,
                                   femtobase::Eta,
                                   femtobase::Phi,
                                   femtotwotrackresonances::ResonanceMass,
                                   femtotwotrackresonances::ResonanceType,
                                   femtotwotrackresonances::PosDauResonanceId,
                                   femtotwotrackresonances::NegDauResonanceId,
                                   femtobase::Theta<femtobase::Eta>,
                                   femtobase::Px<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Py<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                                   femtobase::P<femtobase::Pt, femtobase::Eta>);
using FUResos = FUResos_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUResoMasks_001, "FURESOMASKS", 1,
                                   femtotwotrackresonances::ResonanceMask);
using FUResoMasks = FUResoMasks_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOTWOTRACKRESONANCESDERIVED_H_
