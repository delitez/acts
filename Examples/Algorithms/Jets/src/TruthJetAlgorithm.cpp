// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Jets/TruthJetAlgorithm.hpp"

#include "Acts/Definitions/ParticleData.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/ScopedTimer.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <algorithm>
#include <fstream>
#include <mutex>
#include <ostream>
#include <ranges>
#include <stdexcept>

#include <boost/container/flat_map.hpp>
#include <edm4hep/MCParticle.h>
#include <edm4hep/MCParticleCollection.h>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>

namespace ActsExamples {

TruthJetAlgorithm::TruthJetAlgorithm(const Config& cfg,
                                     Acts::Logging::Level lvl)
    : IAlgorithm("TruthJetAlgorithm", lvl), m_cfg(cfg) {
  if (m_cfg.inputTruthParticles.empty()) {
    throw std::invalid_argument("Input particles collection is not configured");
  }
  m_inputTruthParticles.initialize(m_cfg.inputTruthParticles);
  if (m_cfg.inputTruthParticles.empty()) {
    throw std::invalid_argument("Input particles collection is not configured");
  }
  m_inputEDM4HepParticles.initialize(m_cfg.inputEDM4HepParticles);
  if (m_cfg.inputEDM4HepParticles.empty()) {
    throw std::invalid_argument(
        "Input EDM4Hep particles collection is not configured");
  }
  m_outputJets.initialize(m_cfg.outputJets);
}

namespace {
ActsPlugins::FastJet::JetLabel jetLabelFromHadronType(
    Acts::HadronType hadronType) {
  using enum Acts::HadronType;
  switch (hadronType) {
    case BBbarMeson:
    case BottomMeson:
    case BottomBaryon:
      return ActsPlugins::FastJet::JetLabel::BJet;
    case CCbarMeson:
    case CharmedMeson:
    case CharmedBaryon:
      return ActsPlugins::FastJet::JetLabel::CJet;
    case StrangeMeson:
    case StrangeBaryon:
    case LightMeson:
    case LightBaryon:
      return ActsPlugins::FastJet::JetLabel::LightJet;
    default:
      return ActsPlugins::FastJet::JetLabel::Unknown;
  }
}

}  // namespace

ProcessCode ActsExamples::TruthJetAlgorithm::initialize() {
  if (m_cfg.debugCsvOutput) {
    std::ofstream outfile;
    outfile.open("particles.csv");
    outfile << "event,pt,eta,phi,pdg,label" << std::endl;

    outfile.flush();
    outfile.close();

    outfile.open("jets.csv");
    outfile << "event,pt,eta,phi,label" << std::endl;

    outfile.flush();
    outfile.close();

    outfile.open("hadrons.csv");
    outfile << "event,pt,eta,phi,pdg,label" << std::endl;

    outfile.flush();
    outfile.close();
  }
  return ProcessCode::SUCCESS;
}

ProcessCode ActsExamples::TruthJetAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // Initialize the output container
  std::vector<ActsPlugins::FastJet::TruthJet<TrackContainer>>
      outputJetContainer{};

  Acts::ScopedTimer globalTimer("TruthJetAlgorithm", logger(),
                                Acts::Logging::DEBUG);

  const SimParticleContainer& truthParticlesRaw = m_inputTruthParticles(ctx);
  std::vector<const SimParticle*> truthParticles;
  truthParticles.reserve(truthParticlesRaw.size());
  std::ranges::transform(truthParticlesRaw, std::back_inserter(truthParticles),
                         [](const auto& particle) { return &particle; });
  ACTS_DEBUG("Number of truth particles: " << truthParticles.size());

  const fastjet::JetDefinition defaultJetDefinition = fastjet::JetDefinition(
      fastjet::antikt_algorithm, m_cfg.jetClusteringRadius);

  // Get the 4-momentum information from the simulated truth particles
  // and create fastjet::PseudoJet objects
  std::vector<fastjet::PseudoJet> inputPseudoJets;

  static std::mutex mtxPseudoJets;
  {
    std::ofstream outfile;
    Acts::ScopedTimer timer("Input particle building", logger(),
                            Acts::Logging::DEBUG);

    std::unique_lock lock(mtxPseudoJets, std::defer_lock);
    if (m_cfg.debugCsvOutput) {
      lock.lock();
      outfile.open("particles.csv",
                   std::ios_base::app);  // append instead of overwrite
    }

    for (unsigned int i = 0; i < truthParticles.size(); ++i) {
      const auto* particle = truthParticles.at(i);

      // Check if the generated particle is from the hard scatter, if not, skip
      // Convention is that idx = 0 is the hard scatter
      edm4hep::MutableMCParticle edm4hepParticle;
      EDM4hepUtil::writeParticle(*particle, edm4hepParticle);

      // TO-DO
      //  if (m_cfg.clusterHSParticlesOnly && gp != nullptr &&
      //  HepMC3Util::eventGeneratorIndex(*gp) != 0){
      //    continue;
      //  }

      fastjet::PseudoJet pseudoJet(
          particle->momentum().x(), particle->momentum().y(),
          particle->momentum().z(), particle->energy());

      if (m_cfg.debugCsvOutput) {
        outfile << ctx.eventNumber << "," << pseudoJet.pt() << ","
                << pseudoJet.eta() << "," << pseudoJet.phi() << ","
                << static_cast<int>(particle->pdg()) << ","
                << static_cast<int>(jetLabelFromHadronType(
                       Acts::ParticleIdHelper::hadronType(particle->pdg())));
        outfile << std::endl;
      }

      pseudoJet.set_user_index(i);
      inputPseudoJets.push_back(pseudoJet);
    }

    if (m_cfg.debugCsvOutput) {
      outfile.flush();
      outfile.close();
    }
  }

  ACTS_DEBUG("Number of input pseudo jets: " << inputPseudoJets.size());

  std::vector<fastjet::PseudoJet> jets;
  fastjet::ClusterSequence clusterSeq;
  {
    Acts::ScopedTimer timer("Jet clustering", logger(), Acts::Logging::DEBUG);
    // Run the jet clustering
    clusterSeq =
        fastjet::ClusterSequence(inputPseudoJets, defaultJetDefinition);
    // Get the jets above a certain pt threshold
    jets = sorted_by_pt(clusterSeq.inclusive_jets(m_cfg.jetPtMin));
    // Apply eta range cut if specified
    if (m_cfg.jetEtaRange.first.has_value() ||
        m_cfg.jetEtaRange.second.has_value()) {
      double minEta = m_cfg.jetEtaRange.first.value_or(
          std::numeric_limits<double>::lowest());
      double maxEta =
          m_cfg.jetEtaRange.second.value_or(std::numeric_limits<double>::max());
      std::erase_if(jets, [minEta, maxEta](const auto& jet) {
        return jet.eta() < minEta || jet.eta() > maxEta;
      });
    }
    ACTS_DEBUG("Number of clustered jets: " << jets.size());
  }
  // TO-DO
  // std::vector<std::pair<ActsPlugins::FastJet::JetLabel, std::shared_ptr<const
  // HepMC3::GenParticle>>> hadrons; if(m_cfg.doJetLabeling) {
  //   Acts::ScopedTimer timer("hadron finding", logger(),
  //   Acts::Logging::DEBUG); ACTS_DEBUG("Jet labeling is enabled. Finding
  //   hadrons for jet labeling.");

  //   // const auto& genEvent = *m_inputHepMC3Event(ctx);
  //   const auto& edm4hepParticles = m_inputEDM4HepParticles(ctx);

  //   // A lazy view over the generated particles for hadron finding
  //   auto hadronView =
  //     genEvent.particles() | std::views::filter([this](const auto& particle)
  //     {
  //       if(m_cfg.jetLabelingHSHadronsOnly) {
  //         if (HepMC3Util::eventGeneratorIndex(*particle) != 0) {
  //           return false;
  //         }
  //       }

  //       Acts::PdgParticle pdgId{particle->pdg_id()};
  //       if (!Acts::ParticleIdHelper::isHadron(pdgId)) {
  //         return false;
  //       }

  //       if (particle->status() != HepMC3Util::kDecayedParticleStatus &&
  //           particle->status() != HepMC3Util::kUndecayedParticleStatus) {
  //         return false;
  //       }

  //       // Apply a pt cut on B or C hadrons
  //       auto label = jetLabelFromHadronType(
  //           Acts::ParticleIdHelper::hadronType(pdgId));
  //       using enum ActsPlugins::FastJet::JetLabel;

  //       if(label == BJet || label == CJet) {
  //         if(particle->momentum().pt() < m_cfg.jetLabelingHadronPtMin) {
  //           return false;
  //         }
  //       }
  //       return true;
  //     }) |
  //     std::views::transform([](const auto& particle) {
  //       Acts::PdgParticle pdgId{particle->pdg_id()};
  //       auto type = Acts::ParticleIdHelper::hadronType(pdgId);
  //       auto label = jetLabelFromHadronType(type);
  //       return std::make_pair(label, particle);
  //     }) |
  //     std::views::filter([](const auto& hadron){
  //       return hadron.first > ActsPlugins::FastJet::JetLabel::Unknown;
  //     });

  //     std::ranges::copy(hadronView, std::back_inserter(hadrons));

  //     //Deduplicate hadrons based on their pdg id
  //     std::ranges::sort(hadrons, [](const auto& a, const auto& b){
  //       return a.second->pdg_id() < b.second->pdg_id();
  //     });

  //     auto unique = std::ranges::unique(hadrons);
  //     hadrons.erase(unique.begin(), unique.end());

  //       if (m_cfg.debugCsvOutput) {
  //       static std::mutex mtxHadrons;
  //       std::lock_guard lock(mtxHadrons);
  //       std::ofstream outfile;
  //       outfile.open("hadrons.csv", std::ios_base::app);
  //       for (const auto& hadron : hadrons) {
  //         outfile << ctx.eventNumber << "," << hadron.second->momentum().pt()
  //                 << "," << hadron.second->momentum().eta() << ","
  //                 << hadron.second->momentum().phi() << ","
  //                 << static_cast<int>(hadron.second->pdg_id()) << ","
  //                 << static_cast<int>(hadron.first);
  //         outfile << std::endl;
  //       }
  //     }

  // } // if do jet labeling

  // Jet classification
  // TO-DO
  // auto classifyJet = [&] (const fastjet::PseudoJet& jet) {
  //   auto hadronsInJetView =
  //     hadrons | std::views::filter([&jet, this](const auto& hadron) {
  //       const auto& momentum = hadron.second->momentum();
  //       Acts::Vector3 hadronJetMom{momentum.px(), momentum.py(),
  //       momentum.pz()}; Acts::Vector3 jetMom{jet.px(), jet.py(), jet.pz()};
  //       return Acts::VectorHelpers::deltaR(jetMom, hadronJetMom) <
  //       m_cfg.jetLabelingDeltaR;
  //     }) |
  //     std::views::transform([](const auto& hadron) {
  //       return std::pair{hadron.second,
  //       jetLabelFromHadronType(Acts::ParticleIdHelper::hadronType(
  //           Acts::PdgParticle{hadron.second->pdg_id()}))};
  //     });

  //     std::vector<std::pair<std::shared_ptr<const HepMC3::GenParticle>,
  //     ActsPlugins::FastJet::JetLabel>> hadronsInJet;
  //     std::ranges::copy(hadronsInJetView, std::back_inserter(hadronsInJet));

  //     ACTS_VERBOSE("-> hadrons in jet: " << hadronsInJet.size());
  //     for (const auto& hadron : hadronsInJet) {
  //       ACTS_VERBOSE(
  //           "  - " << hadron.first->pdg_id() << " "
  //                  <<
  //                  Acts::findName(Acts::PdgParticle{hadron.first->pdg_id()}).value_or("UNKNOWN")
  //                  << " label=" << hadron.second);
  //     }

  //     auto maxHadronIt = std::ranges::max_element(hadronsInJet, [](const
  //     auto& a, const auto& b) { return a < b; },
  //       [](const auto& a) {
  //         const auto& [hadron, label] = a;
  //         return label;
  //       });

  //   if (maxHadronIt == hadronsInJet.end()) {
  //     // Now hadronic "jet"
  //     return ActsPlugins::FastJet::JetLabel::Unknown;
  //   }

  //   const auto& [maxHadron, maxHadronLabel] = *maxHadronIt;

  //       ACTS_VERBOSE("-> max hadron type="
  //                <<
  //                Acts::findName(Acts::PdgParticle{maxHadron->pdg_id()}).value_or("UNKNOWN")
  //                << " label=" << maxHadronLabel);

  //   return maxHadronLabel;

  // }; // jet classification

  boost::container::flat_map<ActsPlugins::FastJet::JetLabel, std::size_t>
      jetLabelCounts;

  static std::mutex mtxJets;
  {
    Acts::AveragingScopedTimer timer("Jet classification", logger(),
                                     Acts::Logging::DEBUG);

    std::ofstream outfile;
    std::unique_lock lock(mtxJets, std::defer_lock);
    if (m_cfg.debugCsvOutput) {
      lock.lock();
      outfile.open("jets.csv", std::ios_base::app);
    }

    for (unsigned int i = 0; i < jets.size(); ++i) {
      const auto& jet = jets.at(i);

      // If jet labeling is enabled, classify the jet based on its hadronic
      // content
      ActsPlugins::FastJet::JetLabel jetLabel =
          ActsPlugins::FastJet::JetLabel::Unknown;
      if (m_cfg.doJetLabeling) {
        ACTS_DEBUG("Classifying jet " << i);
        auto sample = timer.sample();
        // TO-DO
        jetLabel = ActsPlugins::FastJet::JetLabel::Unknown;
        // jetLabel = classifyJet(jet);
      }

      if (m_cfg.debugCsvOutput) {
        outfile << ctx.eventNumber << "," << jet.pt() << "," << jet.eta() << ","
                << jet.phi() << "," << static_cast<int>(jetLabel);
        outfile << std::endl;
      }

      // Initialize truth jet for storing in the output container
      Acts::Vector4 jetFourMom{jet.px(), jet.py(), jet.pz(), jet.e()};
      ActsPlugins::FastJet::TruthJet<TrackContainer> truthJet(jetFourMom,
                                                              jetLabel);

      std::vector<fastjet::PseudoJet> jetConstituents = jet.constituents();
      std::vector<int> constituentIndices;
      constituentIndices.reserve(jetConstituents.size());

      for (unsigned int j = 0; j < jetConstituents.size(); ++j) {
        constituentIndices.push_back(jetConstituents[j].user_index());
      }

      truthJet.setConstituentIndices(constituentIndices);
      outputJetContainer.push_back(truthJet);

      jetLabelCounts[jetLabel] += 1;

      ACTS_VERBOSE("-> jet label: " << jetLabel);
      ACTS_VERBOSE("-> jet constituents: ");

      if (logger().doPrint(Acts::Logging::VERBOSE)) {
        for (const auto& constituent : constituentIndices) {
          const auto& particle = truthParticles.at(constituent);
          ACTS_VERBOSE("- " << particle);
        }
      }
    }

    if (m_cfg.debugCsvOutput) {
      outfile.flush();
      outfile.close();
    }
  }

  ACTS_DEBUG("-> jet label counts: ");
  for (const auto& [label, count] : jetLabelCounts) {
    ACTS_DEBUG("  - " << label << ": " << count);
  }

  m_outputJets(ctx, std::move(outputJetContainer));

  return ProcessCode::SUCCESS;
}

ProcessCode ActsExamples::TruthJetAlgorithm::finalize() {
  ACTS_INFO("Finalizing truth jet clustering");
  return ProcessCode::SUCCESS;
}

}  // namespace ActsExamples
