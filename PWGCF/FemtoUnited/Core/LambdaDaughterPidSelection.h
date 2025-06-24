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

/// \file LambdaDaughterPidSelection.h
/// \brief track pid selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_LAMBDADAUGHTERPIDSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_LAMBDADAUGHTERPIDSELECTION_H_

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

#include <cmath>

namespace o2::analysis::femtounited
{
namespace lambdadaughterpidselection
{

// configurables for lambda daughterpid is defined in lambdaselection

/// The different selections this task is capable of doing
enum LambdaDaughterPidSels {

  // charge combination for lambda

  // charge combination for antilambda

  kLambdaDaughterPidSelsMax
};

/// \class LambdaDaughterPidSelectionPos
/// \brief Cut class to contain and execute all cuts applied to pid of vzerodaughters
class LambdaDaughterPidSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackPidMaskType, kLambdaDaughterPidSelsMax>
{
 public:
  LambdaDaughterPidSelection() {}
  virtual ~LambdaDaughterPidSelection() = default;
  template <typename track>
  void applySelections(track const& posDaughter, track const& negDaughter)
  {
    this->reset();

    // apply bits for lambda
    this->setBitmaskForObservable(kNegDaughTpcPion, negDaughter.tpcNSigmaPi());
    this->setBitmaskForObservable(kPosDaughTpcProton, posDaughter.tpcNSigmaPr());

    // apply bits for antilambda
    this->setBitmaskForObservable(kPosDaughTpcPion, posDaughter.tpcNSigmaPi());
    this->setBitmaskForObservable(kNegDaughTpcProton, negDaughter.tpcNSigmaPr());

    this->assembleBismask();
  };

  bool checkDaughterPidLambda()
  {
    return this->getAnySelection(kNegDaughTpcPion) && this->getAnySelection(kPosDaughTpcProton);
  }

  bool checkDaughterPidAntiLambda()
  {
    return this->getAnySelection(kPosDaughTpcPion) && this->getAnySelection(kNegDaughTpcPion);
  }

}; // namespace femtoDream
}; // namespace lambdadaughterpidselection
}; // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_LAMBDADAUGHTERPIDSELECTION_H_
