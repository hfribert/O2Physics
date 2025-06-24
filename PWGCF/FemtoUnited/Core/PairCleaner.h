// Copyright 2019-2022 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file PairCleaner.h
/// \brief pair cleaner class
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_PAIRCLEANER_H_
#define PWGCF_FEMTOUNITED_CORE_PAIRCLEANER_H_

#include "PWGCF/FemtoUnited/Core/DataTypes.h"
#include "PWGCF/FemtoUnited/Core/HistManager.h"
#include "PWGCF/FemtoUnited/Core/Modes.h"

#include "Framework/HistogramRegistry.h"

namespace o2::analysis::femtounited
{
namespace paircleaner
{
/// \class FemtoDreamEventHisto
/// \brief Class for histogramming event properties
// template <femtomodes::Mode mode>
class PairCleaner
{
 public:
  /// Destructor
  virtual ~PairCleaner() = default;
  /// Initializes histograms for the task
  /// \param registry Histogram registry to be passed
  ///
  template <typename T1, typename T2>
  bool isCleanTrackPair(T1 track1, T2 track2)
  {
    return track1.globalIndex() != track2.globalIndex();
  };

 private:
};
}; // namespace paircleaner
}; // namespace o2::analysis::femtounited
#endif // PWGCF_FEMTOUNITED_CORE_PAIRCLEANER_H_
