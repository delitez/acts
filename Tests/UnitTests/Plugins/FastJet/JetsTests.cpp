// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/EventData/AnyTrackProxy.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "ActsPlugins/FastJet/Jets.hpp"

using namespace Acts;
using namespace ActsPlugins;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(FastJetSuite)

BOOST_AUTO_TEST_CASE(TruthParticleOneJet) {
  ActsFatras::Barcode barcode;
  ActsPlugins::FastJet::TruthJet jet(Acts::Vector4(100, 0, 0, 100),
                                     ActsPlugins::FastJet::JetLabel::Unknown);
  jet.setConstituents(std::vector<ActsFatras::Barcode>{barcode});
  BOOST_CHECK_EQUAL(jet.constituents().size(), 1);
  BOOST_CHECK_EQUAL(jet.constituents()[0], barcode);
}

BOOST_AUTO_TEST_CASE(SingleTrackJet) {
  TrackContainer tracks{VectorTrackContainer{}, VectorMultiTrajectory{}};
  auto track = tracks.makeTrack();
  Acts::AnyConstTrackProxy constTrack(track);

  ActsPlugins::FastJet::TrackJet jet(Acts::Vector4(100, 0, 0, 100),
                                     ActsPlugins::FastJet::JetLabel::Unknown);
  std::vector<Acts::AnyConstTrackProxy> constituents{constTrack};
  jet.setConstituents(constituents);

  BOOST_CHECK_EQUAL(jet.constituents().size(), 1);
  BOOST_CHECK(jet.constituents()[0].index() == constTrack.index());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
