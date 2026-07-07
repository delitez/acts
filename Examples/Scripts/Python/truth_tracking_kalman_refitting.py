#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Optional

import acts
import acts.examples
from acts.examples.root import (
    RootTrackStatesWriter,
    RootTrackSummaryWriter,
    RootTrackFitterPerformanceWriter,
)

from truth_tracking_kalman import runTruthTrackingKalman

u = acts.UnitConstants


def runRefittingKf(
    trackingGeometry: acts.TrackingGeometry,
    field: acts.MagneticFieldProvider,
    digiConfigFile: Path,
    outputDir: Path,
    inputHitsPath: Optional[list[Path]] = None,
    multipleScattering: bool = True,
    energyLoss: bool = True,
    reverseFilteringMomThreshold=float("inf"),
    reverseFilteringCovarianceScaling=100.0,
    useJosephFormulation: bool = False,
    s: acts.examples.Sequencer = None,
    args: argparse.Namespace = None,
):
    outputDir.mkdir(parents=True, exist_ok=True)
    s = runTruthTrackingKalman(
        trackingGeometry,
        field,
        digiConfigFile=digiConfigFile,
        outputDir=outputDir,
        inputHitsPath=inputHitsPath,
        reverseFilteringMomThreshold=0 * u.GeV,  # use direct smoothing
        reverseFilteringCovarianceScaling=reverseFilteringCovarianceScaling,
        useJosephFormulation=useJosephFormulation,
        s=s,
        args=args,
    )

    kalmanOptions = {
        "multipleScattering": multipleScattering,
        "energyLoss": energyLoss,
        "reverseFilteringMomThreshold": reverseFilteringMomThreshold,
        "reverseFilteringCovarianceScaling": reverseFilteringCovarianceScaling,
        "freeToBoundCorrection": acts.examples.FreeToBoundCorrection(False),
        "level": acts.logging.INFO,
        "chi2Cut": float("inf"),
        "useJosephFormulation": useJosephFormulation,
    }

    s.addAlgorithm(
        acts.examples.RefittingAlgorithm(
            level=acts.logging.INFO,
            inputTracks="kf_tracks",
            outputTracks="kf_refit_tracks",
            initialVarInflation=6 * [100.0],
            fit=acts.examples.makeKalmanFitterFunction(
                trackingGeometry, field, **kalmanOptions
            ),
            beamSpotConstraint=acts.SquareMatrix2(
                [[(0.0125 * u.mm) ** 2, 0], [0, (55.5 * u.mm) ** 2]]
            ),
        )
    )

    s.addAlgorithm(
        acts.examples.TrackTruthMatcher(
            level=acts.logging.INFO,
            inputTracks="kf_refit_tracks",
            inputParticles="particles_selected",
            inputMeasurementParticlesMap="measurement_particles_map",
            outputTrackParticleMatching="refit_track_particle_matching",
            outputParticleTrackMatching="refit_particle_track_matching",
        )
    )

    # s.addWriter(
    #     RootTrackStatesWriter(
    #         level=acts.logging.INFO,
    #         inputTracks="kf_refit_tracks",
    #         inputParticles="particles_selected",
    #         inputTrackParticleMatching="refit_track_particle_matching",
    #         inputSimHits="simhits",
    #         inputMeasurementSimHitsMap="measurement_simhits_map",
    #         filePath=str(outputDir / "trackstates_kf_refit.root"),
    #     )
    # )

    # s.addWriter(
    #     RootTrackSummaryWriter(
    #         level=acts.logging.INFO,
    #         inputTracks="kf_refit_tracks",
    #         inputParticles="particles_selected",
    #         inputTrackParticleMatching="refit_track_particle_matching",
    #         filePath=str(outputDir / "tracksummary_kf_refit.root"),
    #     )
    # )

    s.addWriter(
        RootTrackFitterPerformanceWriter(
            level=acts.logging.INFO,
            inputTracks="kf_refit_tracks",
            inputParticles="particles_selected",
            inputTrackParticleMatching="refit_track_particle_matching",
            filePath=str(outputDir / "performance_kf_refit.root"),
        )
    )

    return s


def collectEdm4hepInputs(inputPaths: list[Path]) -> list[tuple[Path, Path]]:
    resolvedInputs: list[tuple[Path, Path]] = []

    for inputPath in inputPaths:
        if inputPath.is_dir():
            for edm4hepFile in sorted(inputPath.rglob("edm4hep.root")):
                if edm4hepFile.is_file():
                    relativeParent = edm4hepFile.parent.relative_to(inputPath)
                    outputSuffix = (
                        relativeParent
                        if relativeParent != Path(".")
                        else edm4hepFile.parent.name
                    )
                    resolvedInputs.append((edm4hepFile, outputSuffix))
        else:
            outputSuffix = (
                inputPath.parent.name if inputPath.parent.name else inputPath.stem
            )
            resolvedInputs.append((inputPath, Path(outputSuffix)))

    return resolvedInputs


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run truth tracking Kalman refitting on one or more EDM4hep files"
    )
    parser.add_argument(
        "--edm4hep",
        nargs="+",
        type=Path,
        default=[None],
        help="One or more EDM4hep input files or directories containing edm4hep.root files",
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

    # Get detector
    detector = getOpenDataDetector(odd_dir=geoDir, materialDecorator=oddMaterialDeco)
    trackingGeometry = detector.trackingGeometry()
    field = detector.field

    runRefittingKf(
        trackingGeometry=trackingGeometry,
        field=field,
        digiConfigFile=oddDigiConfig,
        outputDir=args.output,
        # inputHitsPath=args.edm4hep,
        args=args,
    ).run()

    # field = acts.SolenoidBField(
    #     radius=1200 * u.mm,
    #     length=6000 * u.mm,
    #     bMagCenter=3 * u.T,
    #     nCoils=1194,
    # )

    # solenoid = acts.SolenoidBField(
    #     radius=1200 * u.mm,
    #     length=6000 * u.mm,
    #     bMagCenter=3 * u.T,
    #     nCoils=1194,
    # )

    # field = acts.solenoidFieldMap(
    #     rlim=(0, 1200 * u.mm),
    #     zlim=(-5000 * u.mm, 5000 * u.mm),
    #     nbins=(50, 50),
    #     field=solenoid,
    # )
