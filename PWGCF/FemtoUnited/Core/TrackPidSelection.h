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

/// \file TrackPidSelection.h
/// \brief track pid selection
/// \author anton.riedel@tum.de, TU München, anton.riedel@tum.de

#ifndef PWGCF_FEMTOUNITED_CORE_TRACKPIDSELECTION_H_
#define PWGCF_FEMTOUNITED_CORE_TRACKPIDSELECTION_H_

#include "PWGCF/FemtoUnited/Core/BaseSelection.h"
#include "PWGCF/FemtoUnited/Core/DataTypes.h"

#include "Framework/Configurable.h"

#include <cmath>
#include <string>

namespace o2::analysis::femtounited
{
namespace trackpidselection
{


/// The different selections this task is capable of doing
enum TrackPidSels {

  kItsElectron, ///< ITS Electon PID
  kItsPion,     ///< ITS Pion PID
  kItsKaon,     ///< ITS Kaon PID
  kItsProton,   ///< ITS Proton PID
  kItsDeuteron, ///< ITS Deuteron PID
  kItsTriton,   ///< ITS Triton PID
  kItsHelium,   ///< ITS He3 PID

  kTpcElectron, ///< TPC Electon PID
  kTpcPion,     ///< TPC Pion PID
  kTpcKaon,     ///< TPC Kaon PID
  kTpcProton,   ///< TPC Proton PID
  kTpcDeuteron, ///< TPC Deuteron PID
  kTpcTriton,   ///< TPC Triton PID
  kTpcHelium,   ///< TPC He3 PID

  kTofElectron, ///< TOF Electon PID
  kTofPion,     ///< TOF Pion PID
  kTofKaon,     ///< TOF Kaon PID
  kTofProton,   ///< TOF Proton PID
  kTofDeuteron, ///< TOF Deuteron PID
  kTofTriton,   ///< TOF Triton PID
  kTofHelium,   ///< TOF He3 PID

  kTpctofElectron, ///< TPC+TOF Electon PID
  kTpctofPion,     ///< TPC+TOF Pion PID
  kTpctofKaon,     ///< TPC+TOF Kaon PID
  kTpctofProton,   ///< TPC+TOF Proton PID
  kTpctofDeuteron, ///< TPC+TOF Deuteron PID
  kTpctofTriton,   ///< TPC+TOF Triton PID
  kTpctofHelium,   ///< TPC+TOF He3 PID

  kTrackPidSelMax
};

/// \class FemtoDreamTrackCuts
/// \brief Cut class to contain and execute all cuts applied to tracks
class TrackPidSelection : public BaseSelection<float, o2::aod::femtodatatypes::TrackPidMaskType, kTrackPidSelMax>
{
 public:
  TrackPidSelection() {}
  virtual ~TrackPidSelection() = default;
  template <typename track>
  void applySelections(track const& Track)
  {
    this->reset();

    // its pid
    this->evaluateObservable(kItsElectron, Track.itsNSigmaEl());
    this->evaluateObservable(kItsPion, Track.itsNSigmaPi());
    this->evaluateObservable(kItsKaon, Track.itsNSigmaKa());
    this->evaluateObservable(kItsProton, Track.itsNSigmaPr());
    this->evaluateObservable(kItsDeuteron, Track.itsNSigmaDe());
    this->evaluateObservable(kItsTriton, Track.itsNSigmaTr());
    this->evaluateObservable(kItsHelium, Track.itsNSigmaHe());

    // tpc pid
    this->evaluateObservable(kTpcElectron, Track.tpcNSigmaEl());
    this->evaluateObservable(kTpcPion, Track.tpcNSigmaPi());
    this->evaluateObservable(kTpcKaon, Track.tpcNSigmaKa());
    this->evaluateObservable(kTpcProton, Track.tpcNSigmaPr());
    this->evaluateObservable(kTpcDeuteron, Track.tpcNSigmaDe());
    this->evaluateObservable(kTpctofTriton, Track.tpcNSigmaTr());
    this->evaluateObservable(kTpcHelium, Track.tpcNSigmaHe());

    // tof pid
    this->evaluateObservable(kTofElectron, Track.tofNSigmaEl());
    this->evaluateObservable(kTofPion, Track.tofNSigmaPi());
    this->evaluateObservable(kTofKaon, Track.tofNSigmaKa());
    this->evaluateObservable(kTofProton, Track.tofNSigmaPr());
    this->evaluateObservable(kTofDeuteron, Track.tofNSigmaDe());
    this->evaluateObservable(kTofTriton, Track.tofNSigmaTr());
    this->evaluateObservable(kTofHelium, Track.tofNSigmaHe());

    // combined tpc + tof pid
    this->evaluateObservable(kTpctofElectron, std::hypot(Track.tpcNSigmaEl(), Track.tofNSigmaEl()));
    this->evaluateObservable(kTpctofPion, std::hypot(Track.tpcNSigmaPi(), Track.tofNSigmaPi()));
    this->evaluateObservable(kTpctofKaon, std::hypot(Track.tpcNSigmaKa(), Track.tofNSigmaKa()));
    this->evaluateObservable(kTpctofProton, std::hypot(Track.tpcNSigmaPr(), Track.tofNSigmaPr()));
    this->evaluateObservable(kTpctofDeuteron, std::hypot(Track.tpcNSigmaDe(), Track.tofNSigmaDe()));
    this->evaluateObservable(kTpctofTriton, std::hypot(Track.tpcNSigmaTr(), Track.tofNSigmaTr()));
    this->evaluateObservable(kTpctofHelium, std::hypot(Track.tpcNSigmaHe(), Track.tofNSigmaHe()));

    this->assembleBitmask();
  };
}; // namespace femtoDream
}; // namespace trackpidselection
}; // namespace o2::analysis::femtounited

#endif // PWGCF_FEMTOUNITED_CORE_TRACKPIDSELECTION_H_
