from pathlib import Path

import numpy as np
import pandas as pd
import ROOT

ROOT.gROOT.SetBatch(True)  # no canvases pop up during fitting

RootFile = f"/home/jona/phd/key4hep/data/digi/2026-09-24/histograms_NOpt_A_1_B_1_1p2V_0T_singlePlane_20um_hepevt_e-_5.0GeV_cosTheta_8-172deg_G4-EMZ-sDep0-ECut0keV-rCut0.003mm.root"
RootDir = f"VTXBdigi_inner/Layer0/digiHit_residuals"

OutRootFile = Path(__file__).resolve().parent / "cluster_uncertainty_fits.root"


def FetchRootObj(FilePath, Dir, ObjName):
    rFile = ROOT.TFile(FilePath)
    if not rFile:
        raise KeyError("File not found at " + FilePath)

    if Dir:
        rDir = rFile.Get(Dir)
        if not rDir:
            rFile.Close()
            raise KeyError("Directory", Dir, "not found in .root file", FilePath)
    else:
        # No directory given: object lives at the top level of the file
        rDir = rFile

    rObj = rDir.Get(ObjName)
    if not rObj:
        rFile.Close()
        raise KeyError(
            f"Object {ObjName} not found in directory{Dir} in .root file {FilePath}"
        )

    # change ownership of the object to the current process (so it doesn't get deleted when the file is closed)
    rObjClone = rObj.Clone()
    rObjClone.SetDirectory(0)
    rFile.Close()
    return rObjClone


def CalcStdDevTruncated(TH1, Range_NStdDev=2.0):
    Mean = TH1.GetMean()
    StdDev = TH1.GetStdDev()

    RangeMin = Mean - Range_NStdDev * StdDev
    RangeMax = Mean + Range_NStdDev * StdDev

    TH1_local = TH1.Clone()
    TH1_local.SetAxisRange(RangeMin, RangeMax)
    return TH1_local.GetStdDev()


def CalcGausFitSigma(TH1, FitRange_NStdDev=2.0):
    Mean = TH1.GetMean()
    StdDev = TH1.GetStdDev()

    FitMin = Mean - FitRange_NStdDev * StdDev
    FitMax = Mean + FitRange_NStdDev * StdDev

    if TH1.GetEntries() < 2:
        return np.nan
    FitResult = TH1.Fit("gaus", "QSR", "", FitMin, FitMax)
    # Q: quiet, S: return a TFitResult, 0: do not draw the fitted function
    if not FitResult.Get() or not FitResult.IsValid():
        print(f"WARNING: Gaussian fit failed for {TH1.GetName()}")
        return np.nan
    return FitResult.Parameter(2)  # par 0 = norm, 1 = mean, 2 = sigma


def CalcStdDevCentralFraction(TH1, KeepFraction=0.995):
    NBins = TH1.GetNbinsX()
    BinMin, BinMax = 1, NBins
    # BinMin, BinMax = 0, NBins + 1

    Total = TH1.Integral(BinMin, BinMax)
    KeptEntries = Total

    if Total <= 0:
        print(f"WARNING: empty histogram {TH1.GetName()}")
        return np.nan

    TargetEntries = KeepFraction * Total

    while BinMin < BinMax:
        Cost = TH1.GetBinContent(BinMin) + TH1.GetBinContent(BinMax)
        if KeptEntries - Cost <= TargetEntries:
            break

        BinMin += 1
        BinMax -= 1
        KeptEntries -= Cost

    TH1_local = TH1.Clone()
    TH1_local.GetXaxis().SetRange(BinMin, BinMax)
    return TH1_local.GetStdDev()


ClusterLengths = ["1", "2", "3", "4", "5plus"]
Axes = ["u", "v"]

Uncertainties = np.zeros(2 * len(ClusterLengths))
FittedHists = []

for i_ClusterLength, ClusterLength in enumerate(ClusterLengths):
    for i_Axis, Axis in enumerate(Axes):
        RootObj = f"residual_{Axis}_toPrimariesSecondaries_length{ClusterLength}"

        TH1 = FetchRootObj(RootFile, RootDir, RootObj)

        Resolution = TH1.GetStdDev()
        # Resolution = CalcStdDevTruncated(TH1, 2.0)
        # Resolution = CalcGausFitSigma(TH1, 2.0)
        Resolution = CalcStdDevCentralFraction(TH1, 0.95)

        # Uncertainties[i_Axis * len(ClusterLengths) + i_ClusterLength] = TH1.GetStdDev()
        # Uncertainties[i_Axis * len(ClusterLengths) + i_ClusterLength] = (
        #     CalcStdDevTruncated(TH1, 2)
        # )

        Uncertainties[i_Axis * len(ClusterLengths) + i_ClusterLength] = Resolution
        FittedHists.append(TH1)

Uncertainties = Uncertainties * 0.001  # convert from um to mm
print("[", end="")
for i in range(len(Uncertainties) - 1):
    print(f"{Uncertainties[i]:.4f}, ", end="")
print(f"{Uncertainties[-1]:.4f}]")

print("With flipped u/v:")

Uncertainties = np.roll(Uncertainties, len(ClusterLengths))
print("[", end="")
for i in range(len(Uncertainties) - 1):
    print(f"{Uncertainties[i]:.4f}, ", end="")
print(f"{Uncertainties[-1]:.4f}]")

OutFile = ROOT.TFile(str(OutRootFile), "RECREATE")
OutFile.cd()
for Hist in FittedHists:
    Hist.Write()
OutFile.Close()
print(f"Wrote {len(FittedHists)} fitted histograms to {OutRootFile}")
