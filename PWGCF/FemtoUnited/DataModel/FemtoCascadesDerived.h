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

/// \file FemtoCascadesDerived.h
/// \brief cascades tables tables
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoVzerosDerived.h"

namespace o2::aod
{
namespace femtocascades
{
// columns for Vzero
DECLARE_SOA_COLUMN(CascadeMass, cascadeMass, float);                           //! Mass of Vzero
DECLARE_SOA_COLUMN(CascadeMask, cascadeMask, femtodatatypes::CascadeMaskType); //! Bitmask for Vzero selections

// columns for bachelor particles
DECLARE_SOA_INDEX_COLUMN_FULL(Bachelor, bachelor, int, FUTracks, "_Bachelor");                 //!
DECLARE_SOA_COLUMN(BachelorPt, bachelorPt, float);                                             //!
DECLARE_SOA_COLUMN(BachelorEta, bachelorEta, float);                                           //!
DECLARE_SOA_COLUMN(BachelorPhi, bachelorPhi, float);                                           //!
DECLARE_SOA_COLUMN(BachelorMask, bachelorMask, femtodatatypes::CascadeBachelorMaskType);       //!
DECLARE_SOA_COLUMN(BachelorTpcMask, bachelorTpcMask, femtodatatypes::CascadeBachelorMaskType); //!

// columns for daughter vzero
DECLARE_SOA_INDEX_COLUMN_FULL(Vzero, vzero, int, FUVzeros, "_V0"); //!

} // namespace femtocascades

// table for basic cascade information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUCascades_001, "FUCASCADES", 1,
                                   o2::soa::Index<>,
                                   femtobase::CollisionId,
                                   femtobase::Pt,
                                   femtobase::Eta,
                                   femtobase::Phi,
                                   femtocascades::CascadeMass,
                                   femtobase::Theta<femtobase::Eta>,
                                   femtobase::Px<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Py<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                                   femtobase::P<femtobase::Pt, femtobase::Eta>);
using FUCascades = FUCascades_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUCasMasks_001, "FUCASMASK", 1,
                                   femtocascades::CascadeMask);
using FUCasMasks = FUCasMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUCasBacs_001, "FUCASBACS", 1,
                                   femtocascades::BachelorId,
                                   femtocascades::BachelorPt,
                                   femtocascades::BachelorEta,
                                   femtocascades::BachelorPhi,
                                   femtocascades::BachelorMask,
                                   femtocascades::BachelorTpcMask,
                                   femtobase::Theta<femtobase::Eta>,
                                   femtobase::Px<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Py<femtobase::Pt, femtobase::Eta>,
                                   femtobase::Pz<femtobase::Pt, femtobase::Eta>,
                                   femtobase::P<femtobase::Pt, femtobase::Eta>)
using FUCasBacs = FUCasBacs_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUCasVzeros_001, "FUCASVZEROS", 1,
                                   femtocascades::VzeroId,
                                   femtobase::Pt,
                                   femtobase::Eta,
                                   femtobase::Phi,
                                   femtovzeros::VzeroMass,
                                   femtovzeros::PosDauId,
                                   femtovzeros::PosDauPt,
                                   femtovzeros::PosDauEta,
                                   femtovzeros::PosDauPhi,
                                   femtovzeros::PosDauTrackMask,
                                   femtovzeros::PosDauTPCMask,
                                   femtovzeros::NegDauId,
                                   femtovzeros::NegDauPt,
                                   femtovzeros::NegDauEta,
                                   femtovzeros::NegDauPhi,
                                   femtovzeros::NegDauTrackMask,
                                   femtovzeros::NegDauTPCMask);
using FUCasVzeros = FUCasVzeros_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOCASCADESDERIVED_H_
