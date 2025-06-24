
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

/// \file FemtoBaseDerived.h
/// \brief base tables
/// \author Anton Riedel, TU München, anton.riedel@cern.ch

#ifndef PWGCF_FEMTOUNITED_DATAMODEL_FEMTOBASEDERIVED_H_
#define PWGCF_FEMTOUNITED_DATAMODEL_FEMTOBASEDERIVED_H_

#include "Framework/ASoA.h"
#include "Framework/Expressions.h"

#include <cmath>

namespace o2::aod
{

namespace femtobase
{
// Basic information need for all types of particles
DECLARE_SOA_INDEX_COLUMN(Collision, collision); //! collision index
DECLARE_SOA_COLUMN(Pt, pt, float);              //! pt
DECLARE_SOA_COLUMN(Eta, eta, float);            //! eta
DECLARE_SOA_COLUMN(Phi, phi, float);            //! phi

// dynamic columns
// use fabs for pt so it can also be used with signed pt
DECLARE_SOA_DYNAMIC_COLUMN(Theta, theta, //! Theta
                           [](float eta) -> float {
                             return 2.f * std::atan(std::exp(-eta));
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Px, px, //! Momentum in x
                           [](float pt, float phi) -> float {
                             return std::fabs(pt) * std::sin(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Py, py, //! Momentum in y
                           [](float pt, float phi) -> float {
                             return std::fabs(pt) * std::cos(phi);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(Pz, pz, //! Momentum in z
                           [](float pt, float eta) -> float {
                             return std::fabs(pt) * std::sinh(eta);
                           });
DECLARE_SOA_DYNAMIC_COLUMN(P, p, //! Momentum
                           [](float pt, float eta) -> float {
                             return std::fabs(pt) * std::cosh(eta);
                           });
} // namespace femtobase
} // namespace o2::aod
#endif // PWGCF_FEMTOUNITED_DATAMODEL_FEMTOBASEDERIVED_H_
