// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TrackJets/TruthJetAlgorithm.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <ostream>
#include <stdexcept>

namespace ActsExamples{

TruthJetAlgorithm::TruthJetAlgorithm(const Config& cfg,
                                                 Acts::Logging::Level lvl)
    : IAlgorithm("TruthJetAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputTruthParticles.empty()) {
    throw std::invalid_argument("Input particles collection is not configured");
  }
  m_inputTruthParticles.initialize(m_cfg.inputTruthParticles);
  m_outputJets.initialize(m_cfg.outputJets);
  
}

ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {

  const auto& truthParticles = m_inputTruthParticles(ctx);

  ACTS_DEBUG("Number of truth particles: " << truthParticles.size());

  // Get the 4-momentum information from the simulated truth particles
  // and create fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;

  int particleIndex = 0;
  for (const auto& particle : truthParticles) {

    fastjet::PseudoJet pseudoJet( particle.momentum().x(),
                                  particle.momentum().y(),
                                  particle.momentum().z(),
                                  particle.energy() );

    pseudoJet.set_user_index(particleIndex);
    inputPseudoJets.push_back(pseudoJet);
    particleIndex++;

    }
    ACTS_DEBUG("Number of input pseudo jets: " << inputPseudoJets.size());
    // Run the jet clustering
    fastjet::ClusterSequence clusterSeq(inputPseudoJets, DefaultJetDefinition);
    // Get the jets above pT 20 GeV
    std::vector<fastjet::PseudoJet> jets = sorted_by_pt(clusterSeq.inclusive_jets(10));

    if(m_writeJetRootFile) {
          // Write the jets to a ROOT file
          std::unique_ptr<TFile> file( TFile::Open("jetKinematics.root", "RECREATE") );
            if (!file || file->IsZombie()) { throw std::runtime_error("Failed to open file.root"); }
          auto tree = std::make_unique<TTree>("tree", "Jet Tree");

          std::vector<float> jetPx;
          std::vector<float> jetPy;
          std::vector<float> jetPz;
          std::vector<float> jetPt;
          std::vector<float> jetEta;
          std::vector<float> jetPhi;
          std::vector<float> jetE;
          std::vector<float> jetMass;

          tree->Branch("jetPx", &jetPx);
          tree->Branch("jetPy", &jetPy);
          tree->Branch("jetPz", &jetPz);
          tree->Branch("jetPt", &jetPt);
          tree->Branch("jetEta", &jetEta);
          tree->Branch("jetPhi", &jetPhi);
          tree->Branch("jetE", &jetE);
          tree->Branch("jetMass", &jetMass);
          

          for(int i = 0; i < jets.size(); i++) {
            ACTS_DEBUG("DEBUG: Jet momentum: " << jets[i].px() << ", "
                      << jets[i].py() << ", "
                      << jets[i].pz() << ", "
                      << jets[i].E());
            jetPt.push_back(jets[i].perp());
            jetPx.push_back(jets[i].px());
            jetPy.push_back(jets[i].py());
            jetPz.push_back(jets[i].pz());
            jetE.push_back(jets[i].E());
            jetEta.push_back(jets[i].eta());
            jetPhi.push_back(jets[i].phi());
            jetMass.push_back(jets[i].m());
          }

          tree->Fill();
          tree->Write();
        }
    // Store the jets in the output data handle
    m_outputJets(ctx, std::move(jets));

  return ProcessCode::SUCCESS;
}

ProcessCode ActsExamples::TruthJetAlgorithm::finalize() const {
	ACTS_INFO(
      "Finalizing truth jet clustering");
	return ProcessCode::SUCCESS;
}


}; // namespace ActsExamples
