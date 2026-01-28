#!/usr/bin/env python3

from pathlib import Path
from typing import Optional
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
import sys
import time

import acts
import acts.examples
from acts.examples.edm4hep import EDM4hepParticleOutputConverter, PodioWriter, PodioReader, EDM4hepSimInputConverter
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

actsDir = Path(__file__).parent.parent.parent.parent

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants

from acts.examples.simulation import (
    addParticleGun,
    ParticleConfig,
    EtaConfig,
    PhiConfig,
    MomentumConfig,
    addFatras,
    addPythia8,
    TruthJetConfig,
    addTruthJetAlg,
    addDigitization,
    ParticleSelectorConfig,
    addSimParticleSelection,
    addDigiParticleSelection,
    addGenParticleSelection,
)
from acts.examples.reconstruction import (
    addSeeding,
    SeedingAlgorithm,
    addKalmanTracks,
    addVertexFitting,
    VertexFinder,
)

from acts.examples.root import (
    RootTrackStatesWriter,
    RootTrackSummaryWriter,
    RootTrackFitterPerformanceWriter,
)

outputDir =  Path.cwd() / "test_output_edm4hep_2"
geoDir = getOpenDataDetectorDirectory()

oddMaterialMap = (geoDir / "data/odd-material-maps.root")
oddDigiConfig = (actsDir / "Examples/Configs/odd-digi-smearing-config.json")
oddSeedingSel = actsDir / "Examples/Configs/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

detector = getOpenDataDetector(odd_dir=geoDir, materialDecorator=oddMaterialDeco)
trackingGeometry = detector.trackingGeometry()
decorators = detector.contextDecorators()
field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(
            events=100,
            numThreads=-1,
            logLevel=acts.logging.INFO,
            outputDir=str(outputDir),
            trackFpes=False,
    )

import acts.examples.edm4hep
from acts.examples.edm4hep import PodioReader
# from DDSim.DD4hepSimulation import DD4hepSimulation

# addPythia8(
#         s,
#         nhard=100,
#         npileup=0,
#         # hardProcess=[
#         #     "HardQCD:all = off",
#         #     "HardQCD:gg2bbbar = on",
#         #     "HardQCD:qqbar2bbbar = on",
#         #     "PartonLevel:ISR = off",
#         #     "PartonLevel:FSR = off",
#         #     "PartonLevel:MPI = off",
#         #     "HadronLevel:all = on",
#         #     "511:mayDecay = off",
#         #     "521:mayDecay = off",
#         #     "531:mayDecay = off",
#         #     "541:mayDecay = off",
#         #     "5122:mayDecay = off",
#         # ],
#         # hardProcess=[
#         #     "Top:qqbar2ttbar=on",
#         #     "HadronLevel:Decay=off",
#         #     "StringFlav:probQQtoQ = 0.0",
#         # ],
#         # hardProcess=["WeakSingleBoson:ffbar2gmZ = on","SoftQCD:nonDiffractive=on","HardQCD:hardbbbar=on","HadronLevel:all = on"],
#         # hardProcess=["Top:qqbar2ttbar=on", "HadronLevel:Decay = off","StringFlav:probQQtoQ = 0.0","ParticleDecays:limitTau0=off"],
#         # hardProcess=["WeakSingleBoson:ffbar2gmZ = on","23:onMode = off", "23:onIfAny = 5","ParticleDecays:limitTau0=off"],
#         hardProcess=["Top:qqbar2ttbar=on"],
#         vtxGen=acts.examples.GaussianVertexGenerator(
#             mean=acts.Vector4(0, 0, 0, 0),
#             stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
#         ),
#         rnd=rnd,
#         outputDirRoot=outputDir,
#         outputDirCsv=outputDir,
#         writeHepMC3=outputDir / "hepmc3_particles.hepmc",
#     )


# odd_xml_file = str(
#             getOpenDataDetectorDirectory() / "xml" / "OpenDataDetector.xml"
#         )

# input = odd_xml_file
# print("Generating input file for test_odd_edm4hep:", input)

# ddsim = DD4hepSimulation()
# ddsim.random.seed = 37
# ddsim.compactFile = [input]
# ddsim.enableGun = False
# ddsim.inputFile = str("/Users/delitez/atlas/acts-cvmfs/test_output_2/hepmc3_particles.hepmc")
# #ddsim.inputFile = outputDir / "hepmc3_particles.hepmc"
# ddsim.outputFile = outputDir / "ddsim.edm4hep.root"
# ddsim.outputConfig.forceEDM4HEP = True
# ddsim.run()

s.addReader(
        PodioReader(
            level=acts.logging.INFO,
            inputPath=str(outputDir / "test.edm4hep.root"),
            outputFrame="events",
            category="events",
        )
    )

edm4hepReader = acts.examples.edm4hep.EDM4hepSimInputConverter(
        inputFrame="events",
        inputSimHits=[
            "PixelBarrelReadout",
            "PixelEndcapReadout",
            "ShortStripBarrelReadout",
            "ShortStripEndcapReadout",
            "LongStripBarrelReadout",
            "LongStripEndcapReadout",
        ],
        outputParticlesGenerator="particles_generated",
        outputParticlesSimulation="particles_simulated",
        outputSimHits="simhits",
        outputSimVertices="vertices_truth",
        dd4hepDetector=detector,
        trackingGeometry=trackingGeometry,
        sortSimHitsInTime=False,
        particleRMax=1080 * u.mm,
        particleZ=(-3030 * u.mm, 3030 * u.mm),
        particlePtMin=150 * u.MeV,
        level=acts.logging.INFO,
    )
s.addAlgorithm(edm4hepReader)

s.addWhiteboardAlias("particles", edm4hepReader.config.outputParticlesSimulation)

addSimParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0, 24 * u.mm),
            absZ=(0.0, 1.0 * u.m),
            eta=(-3.0, 3.0),
            removeNeutral=True,
        ),
    )

addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=oddDigiConfig,
        rnd=rnd,
        logLevel=acts.logging.ERROR,
    )

addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_generated",
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        particleHypothesis=acts.ParticleHypothesis.muon,
    )

reverseFilteringMomThreshold = 0 * u.GeV

addKalmanTracks(
        s,
        trackingGeometry,
        field,
        reverseFilteringMomThreshold,
        logLevel=acts.logging.FATAL,
    )

s.addAlgorithm(
        acts.examples.TrackSelectorAlgorithm(
            level=acts.logging.INFO,
            inputTracks="tracks",
            outputTracks="selected-tracks",
            selectorConfig=acts.TrackSelector.Config(
                minMeasurements=7,
            ),
        )
    )
s.addWhiteboardAlias("tracks", "selected-tracks")

addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.AMVF,
        outputDirRoot=outputDir,
        logLevel=acts.logging.FATAL,
    )

addTruthJetAlg(
        s,
        TruthJetConfig(
            inputTruthParticles="particles_generated",
            outputJets="output_jets",
            jetPtMin=10 * u.GeV,
        ),
        loglevel=acts.logging.INFO,
    )

s.addWriter(
        acts.examples.root.RootTrackSummaryWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_generated",
            inputJets="output_jets",
            writeJets=True,
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "tracksummary_kf.root"),
        )
    )

s.run()


