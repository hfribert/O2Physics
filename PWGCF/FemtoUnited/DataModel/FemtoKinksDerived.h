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

/// \file FemtoKinksDerived.h
/// \brief kink tables
/// \author Anton Riedel, TU München, anton.riedel@cern.ch
/// \author Henrik Fribert, TU München, henrik.fribert@cern.ch

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOKINKSDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOKINKSDERIVED_H_

#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoBaseDerived.h"
#include "PWGCF/FemtoUnited/DataModel/FemtoTracksDerived.h"

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"

namespace o2::aod
{
namespace femtokinks
{
// columns for bit masks
DECLARE_SOA_COLUMN(MaskKink, maskKink, femtodatatypes::KinkMaskType); //! Bitmask for kink selections

// columns for debug information
// DECLARE_SOA_COLUMN(MassK0short, massK0short, float); //! Mass of Lambda
// DECLARE_SOA_COLUMN(MassAnti, massAnti, float);       //! Mass of Lambda
DECLARE_SOA_COLUMN(CosPa, cosPa, float);             //! Sigma charged daughter DCA at decay vertex
DECLARE_SOA_COLUMN(DauDCA, dauDca, float);           //! Sigma charged daughter DCA at decay vertex
DECLARE_SOA_COLUMN(TransRadius, transRadius, float); //! Sigma transverse radius
DECLARE_SOA_COLUMN(DecayVtxX, decayVtxX, float);     //! x coordinate of Sigma decay vertex
DECLARE_SOA_COLUMN(DecayVtxY, decayVtxY, float);     //! y coordinate of Sigma decay vertex
DECLARE_SOA_COLUMN(DecayVtxZ, decayVtxZ, float);     //! z coordinate of Sigma decay vertex

// id column for charged daughter track
DECLARE_SOA_INDEX_COLUMN_FULL(ChaDau, chaDau, int32_t, FUTracks, "_ChaDau"); //!
} // namespace femtokinks

// table for basic sigma minus information
DECLARE_SOA_TABLE_STAGED_VERSIONED(FUSigmas_001, "FUSIGMAS", 1,
    o2::soa::Index<>,
    femtobase::stored::CollisionId, // use sign to differentiate between sigma minus (-1) and anti sigma minus (+1)
    femtobase::stored::SignedPt,
    femtobase::stored::Eta,
    femtobase::stored::Phi,
    femtobase::stored::Mass,
    femtokinks::ChaDauId,
    femtobase::dynamic::Sign<femtobase::stored::SignedPt>,
    femtobase::dynamic::Pt<femtobase::stored::SignedPt>,
    femtobase::dynamic::P<femtobase::stored::SignedPt, femtobase::stored::Eta>,
    femtobase::dynamic::Px<femtobase::stored::SignedPt, femtobase::stored::Eta>,
    femtobase::dynamic::Py<femtobase::stored::SignedPt, femtobase::stored::Eta>,
    femtobase::dynamic::Pz<femtobase::stored::SignedPt, femtobase::stored::Eta>,
    femtobase::dynamic::Theta<femtobase::stored::Eta>);
using FUSigmas = FUSigmas_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUSigmaMasks_001, "FUSIGMAMASKS", 1,
    femtokinks::MaskKink);
using FUSigmaMasks = FUSigmaMasks_001;

DECLARE_SOA_TABLE_STAGED_VERSIONED(FUSigmaExtras_001, "FUSIGMAEXTRAS", 1,
    femtokinks::CosPa,
    femtokinks::DauDCA,
    femtokinks::DecayVtxX,
    femtokinks::DecayVtxY,
    femtokinks::DecayVtxZ,
    femtokinks::TransRadius);
    // femtokinks::MassAnti,
    // femtokinks::MassK0short

using FUSigmaExtras = FUSigmaExtras_001;

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOKINKSDERIVED_H_