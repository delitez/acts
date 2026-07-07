#!/usr/bin/env python3

from pathlib import Path
from typing import Optional, Union
import argparse

from acts.examples.odd import getOpenDataDetector
import acts
import acts.examples
from acts.examples.edm4hep import (
    EDM4hepParticleOutputConverter,
    PodioWriter,
    PodioReader,
    EDM4hepSimInputConverter,
)

u = acts.UnitConstants


def runTruthTrackingKalman(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    inputParticlePath: Optional[Path] = None,
    inputHitsPath: Optional[Union[Path, list[Path]]] = None,
    decorators=[],
    generatedParticleType: acts.PdgParticle = acts.PdgParticle.eMuon,
    reverseFilteringMomThreshold=0 * u.GeV,
    reverseFilteringCovarianceScaling=100.0,
    numParticles=1,
    linkForward: bool = False,
    useJosephFormulation: bool = False,
    s: acts.examples.Sequencer = None,
    args: argparse.Namespace = None,
):
    from acts.examples.simulation import (
        addParticleGun,
        ParticleConfig,
        EtaConfig,
        PhiConfig,
        MomentumConfig,
        addGenParticleSelection,
        addFatras,
        addPythia8,
        addGeant4,
        addDigitization,
        addSimParticleSelection,
        ParticleSelectorConfig,
        addDigiParticleSelection,
    )

    from acts.examples.root import (
        RootParticleReader,
        RootSimHitReader,
        RootTrackStatesWriter,
        RootTrackSummaryWriter,
        RootTrackFitterPerformanceWriter,
    )

    from acts.examples.reconstruction import (
        addSeeding,
        SeedingAlgorithm,
        TrackSmearingSigmas,
        addKalmanTracks,
        addCKFTracks,
        TrackSelectorConfig,
        CkfConfig,
    )

    if s is None:
        s = acts.examples.Sequencer(
            events=200,
            numThreads=-1,
            logLevel=acts.logging.INFO,
        )

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    for d in decorators:
        s.addContextDecorator(d)

    rnd = acts.examples.RandomNumbers(seed=42)
    outputDir = Path(outputDir)

    logger = acts.getDefaultLogger("Truth tracking example", acts.logging.INFO)

    if args.edm4hep:
        print(
            "Running truth tracking Kalman refitting on EDM4hep input files:",
            args.edm4hep,
        )
        ## Read both particles and hits from EDM4hep input
        if inputHitsPath is not None:
            inputHitsFiles = (
                [inputHitsPath] if isinstance(inputHitsPath, Path) else inputHitsPath
            )
            detector = getOpenDataDetector()
            # reader takes a vector of paths, so we can pass the list directly.
            s.addReader(
                PodioReader(
                    level=acts.logging.DEBUG,
                    inputPath=[str(path) for path in inputHitsFiles],
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

            s.addWhiteboardAlias(
                "particles", edm4hepReader.config.outputParticlesSimulation
            )

            addSimParticleSelection(
                s,
                ParticleSelectorConfig(
                    rho=(0.0, 24 * u.mm),
                    absZ=(0.0, 1.0 * u.m),
                    eta=(-3.0, 3.0),
                    removeNeutral=True,
                ),
            )
    if args.ttbar:
        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"],
            npileup=200,
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
            ),
            rnd=rnd,
        )
        addGenParticleSelection(
            s,
            ParticleSelectorConfig(
                rho=(0.0, 24 * u.mm),
                absZ=(0.0, 1.0 * u.m),
                eta=(-3.0, 3.0),
                pt=(150 * u.MeV, None),
            ),
        )
    if args.particleGun:
        addParticleGun(
            s,
            ParticleConfig(
                num=numParticles, pdg=generatedParticleType, randomizeCharge=True
            ),
            EtaConfig(-3.0, 3.0, uniform=True),
            MomentumConfig(1.0 * u.GeV, 100.0 * u.GeV, transverse=True),
            PhiConfig(0.0, 360.0 * u.degree),
            vtxGen=acts.examples.GaussianVertexGenerator(
                mean=acts.Vector4(0, 0, 0, 0),
                stddev=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=1,
            rnd=rnd,
        )
    # else:
    #     logger.info("Reading particles from {}", inputParticlePath.resolve())
    #     assert inputParticlePath.exists()
    #     s.addReader(
    #             RootParticleReader(
    #                 level=acts.logging.INFO,
    #                 filePath=str(inputParticlePath.resolve()),
    #                 outputParticles="particles_generated",
    #             )
    #         )
    #     s.addWhiteboardAlias("particles", "particles_generated")

    if args.fullsim:
        if s.config.numThreads != 1:
            raise ValueError("Geant 4 simulation does not support multi-threading")
        detector = getOpenDataDetector()
        addGeant4(
            s,
            detector,
            trackingGeometry,
            field,
            rnd=rnd,
            killVolume=trackingGeometry.highestTrackingVolume,
            killAfterTime=25 * u.ns,
            killSecondaries=True,
            logLevel=acts.logging.INFO,
        )
    if args.fatras:
        addFatras(
            s,
            trackingGeometry,
            field,
            rnd=rnd,
            enableInteractions=True,
        )
    # else:
    #     logger.info("Reading hits from {}", inputHitsPath.resolve())
    #     assert inputHitsPath.exists()
    #     s.addReader(
    #             RootSimHitReader(
    #                 level=acts.logging.INFO,
    #                 filePath=str(inputHitsPath.resolve()),
    #                 outputSimHits="simhits",
    #             )
    #         )
    #     s.addWhiteboardAlias("particles_simulated_selected", "particles_generated")

    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        rnd=rnd,
    )

    addDigiParticleSelection(
        s,
        ParticleSelectorConfig(
            pt=(0.9 * u.GeV, None),
            measurements=(7, None),
            removeNeutral=True,
            removeSecondaries=True,
        ),
    )

    addSeeding(
        s,
        trackingGeometry,
        field,
        rnd=rnd,
        inputParticles="particles_generated",
        seedingAlgorithm=SeedingAlgorithm.TruthSmeared,
        trackSmearingSigmas=TrackSmearingSigmas(
            # zero everything so the KF has a chance to find the measurements
            loc0=0,
            loc0PtA=0,
            loc0PtB=0,
            loc1=0,
            loc1PtA=0,
            loc1PtB=0,
            time=0,
            phi=0,
            theta=0,
            ptRel=0,
        ),
        particleHypothesis=acts.ParticleHypothesis.muon,
        initialSigmas=[
            1 * u.mm,
            1 * u.mm,
            1 * u.degree,
            1 * u.degree,
            0 / u.GeV,
            1 * u.ns,
        ],
        initialSigmaQoverPt=0.1 / u.GeV,
        initialSigmaPtRel=0.1,
        initialVarInflation=[1e0, 1e0, 1e0, 1e0, 1e0, 1e0],
        geoSelectionConfigFile=srcdir / "Examples/Configs/odd-seeding-config.json",
    )

    addKalmanTracks(
        s,
        trackingGeometry,
        field,
        reverseFilteringMomThreshold,
        reverseFilteringCovarianceScaling,
        linkForward=linkForward,
        useJosephFormulation=useJosephFormulation,
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

    # s.addWriter(
    #     RootTrackStatesWriter(
    #         level=acts.logging.INFO,
    #         inputTracks="tracks",
    #         inputParticles="particles_selected",
    #         inputTrackParticleMatching="track_particle_matching",
    #         inputSimHits="simhits",
    #         inputMeasurementSimHitsMap="measurement_simhits_map",
    #         filePath=str(outputDir / "trackstates_kf.root"),
    #     )
    # )

    # s.addWriter(
    #     RootTrackSummaryWriter(
    #         level=acts.logging.INFO,
    #         inputTracks="tracks",
    #         inputParticles="particles_selected",
    #         inputTrackParticleMatching="track_particle_matching",
    #         filePath=str(outputDir / "tracksummary_kf.root"),
    #     )
    # )

    s.addWriter(
        RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="track_particle_matching",
            filePath=str(outputDir / "performance_kf.root"),
        )
    )

    return s


if "__main__" == __name__:
    import argparse

    parser = argparse.ArgumentParser(
        description="Run the truth-tracking Kalman example on one or more EDM4hep files"
    )
    parser.add_argument(
        "--edm4hep",
        nargs="+",
        type=Path,
        default=[Path("/Users/delitez/atlas/ttbar_edm4hep/0/edm4hep.root")],
        help="One or more EDM4hep input files",
    )
    parser.add_argument(
        "--fullsim",
        action="store_true",
        help="Run with Geant4 simulation instead of EDM4hep input",
    )
    parser.add_argument(
        "--fatras",
        action="store_true",
        help="Run with Fatras simulation instead of EDM4hep input",
    )
    parser.add_argument(
        "--ttbar",
        action="store_true",
        help="Run with Pythia8 ttbar events instead of EDM4hep input",
    )
    parser.add_argument(
        "--particleGun",
        action="store_true",
        help="Run with particle gun instead of EDM4hep input",
    )
    parser.add_argument(
        "--output",
        "-o",
        type=Path,
        default=Path.cwd(),
        help="Output directory",
    )

    args = parser.parse_args()

    srcdir = Path(__file__).resolve().parent.parent.parent.parent

    from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory

    # Get detector and field
    geoDir = getOpenDataDetectorDirectory()

    # Load material map
    oddMaterialMap = geoDir / "data/odd-material-maps.root"
    oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"

    oddSeedingSel = geoDir / "config/odd-seeding-config.json"
    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

    oddSeedingSel = geoDir / "config/odd-seeding-config.json"
    oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap)

    # Get detector
    detector = getOpenDataDetector(odd_dir=geoDir, materialDecorator=oddMaterialDeco)
    trackingGeometry = detector.trackingGeometry()
    field = detector.field

    runTruthTrackingKalman(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=oddDigiConfig,
        outputDir=args.output,
        inputHitsPath=args.edm4hep,
        args=args,
    ).run()
