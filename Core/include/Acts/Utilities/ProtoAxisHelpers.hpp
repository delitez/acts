// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/ProtoAxis.hpp"


#include <array>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Acts {

inline std::size_t binsFromProtoAxis(DirectedProtoAxis& axis) {
  return axis.getAxis().getNBins();
}

inline std::size_t binsFromProtoAxes(const std::vector<DirectedProtoAxis>& axes) {
  std::size_t nBins = 0;
  for (std::size_t i = 0; i < axes.size(); i++)
  {
    nBins += axes[i].getAxis().getNBins();
  }
  return nBins;
}

inline std::size_t binFromProtoAxis(DirectedProtoAxis& axis, const Vector2& lp) {
  BinningData bd(axis);
  return bd.searchLocal(lp);
}

inline std::array<std::size_t, 3> binTripleFromProtoAxes(
    const std::vector<DirectedProtoAxis>& axes, const Vector3& gp) {
    // const Transform3& invTranform3 = Transform3::Identity(); 
    // const Vector3& bPosition = gp * invTranform3;  
    const Vector3& bPosition = gp; // @todo apply inverse transform?
    std::array<std::size_t, 3> bTriple = {0, 0, 0};
    if (axes.size() > 0){
      BinningData bd0(axes[0]);
      bTriple[0] = bd0.searchGlobal(bPosition);}
    if (axes.size() > 1){
      BinningData bd1(axes[1]);
      bTriple[1] = bd1.searchGlobal(bPosition);}
    if (axes.size() > 2){
      BinningData bd2(axes[2]);
      bTriple[2] = bd2.searchGlobal(bPosition);}
    return bTriple;
}

} // namespace Acts