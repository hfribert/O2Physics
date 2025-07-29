// Copyright 2019-2024 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file dataTypes.h
/// \brief datatypes for bitmasks
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_DATATYPES_H_
#define PWGCF_FEMTOUNITED_CORE_DATATYPES_H_

#include <cstdint>

namespace o2::aod
{
namespace femtodatatypes
{
// Note: Length of the bitmask is the limit of how many selections can be configured

// datatypes for tracks
using TrackMaskType = uint64_t;
using TrackType = uint16_t;

// datatypes for v0s
using V0MaskType = uint32_t;
using V0Type = uint8_t;

// datatypes for kinks
using KinkMaskType = uint32_t;
using MaskTpye = uint16_t;

// datatypes for two track resonances
using TwoTrackResonanceMaskType = uint64_t;
// two track resonance types
using TwoTrackResonanceType = uint8_t;

// datatypes for cascades
using CascadeMaskType = uint32_t;
using CascadeType = uint8_t;

} // namespace femtodatatypes

} // namespace o2::aod

#endif // PWGCF_FEMTOUNITED_CORE_DATATYPES_H_
