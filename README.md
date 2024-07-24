LeukoCt
================

- [Installation Guide](#installation-guide)
  - [Content](#content)
  - [Software & Hardware Requirements](#software--hardware-requirements)
- [Operation & DEMO](#operation--demo)
  - [`STEP02_to_07_Batch_Pipeline.m`](#step02_to_07_batch_pipelinem)
  - [`WBCConcStabilityAnalysis.m`](#wbcconcstabilityanalysism)
- [Citation](#citation)

**LeukoCt** measures the circulating leukocyte concentration in videos
of the human oral mucosal microvasculature captured by a miniaturized
back-illumination microscope (mOBM). It allows the user to define a
region of interest (ROI), and subsequently reports the leukocyte count,
blood flow velocity and volume, the leukocyte concentration, and whether
the concentration over time is sufficiently stable to be a robust
predictor of the clinical WBCC value.

## Installation Guide

### Content

**LeukoCt** includes 4 subfolders:

1.  `Scripts`: includes .m files to be run by user
2.  `Dependencies`: includes custom written functions, and packages from
    [MATLAB File
    Exchange](https://www.mathworks.com/matlabcentral/fileexchange/)
    used by the software  
3.  `Data_In`: includes sample mOBM videos used for DEMO. Due to
    GitHub’s file size limit, please visit [this
    link](https://www.dropbox.com/scl/fo/0l28s2k0mpbvak1hvxx3n/ADwg2e85o6eFZscNxLrCWN4?rlkey=slz2x244q0mcr5qjgq6lc39r7&dl=0)
    to download the files.  
4.  `Data_Out`: includes the output of DEMO. Due to GitHub’s file size
    limit, please visit [this
    link](https://www.dropbox.com/scl/fo/0l28s2k0mpbvak1hvxx3n/ADwg2e85o6eFZscNxLrCWN4?rlkey=slz2x244q0mcr5qjgq6lc39r7&dl=0)
    to download the files.

After downloading the Zip file, unzip and move the folder into the
MATLAB working directory. Add folder and subfolders to MATLAB path.

### Software & Hardware Requirements

**LeukoCt** was written and tested on MATLAB_R2023b (Version 23.2), and
requires the following toolboxes of the same version:

- Signal Processing Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox
- Parallel Computing Toolbox
- Computer Vision Toolbox

Custom written functions, and packages from [MATLAB File
Exchange](https://www.mathworks.com/matlabcentral/fileexchange/) are
included in the “Dependencies” subfolder. Please see Citation section
below.

A 2023 Apple M2 Max MacBook Pro (Model MPHG3LL/A; 32 GB Memory) was used
to run the software. With parallel pool (with Threads; 12 workers)
turned on, a complete run of the DEMO is expected to take 30 minutes.

## Operation & DEMO

The input of **LeukoCt** are raw videos captured by mOBM. The videos are
processed by scripts in the following order:

- `STEP00_Batch_FastStackReg_XYZIntensityEq.m`
- `STEP01_VesselROI.m`
- `STEP02_SCS_VesselSz.m`
- `STEP03_Velocity.m`
- `STEP04_Velocity_Volume.m`
- `STEP05_WBCCount.m`
- `STEP06_WBCConc.m`
- `STEP07_Report.m`
- `WBCConcStabilityAnalysis.m` <br> <br> \####
  `STEP00_Batch_FastStackReg_XYZIntensityEq` This script registers all
  videos in the folder `LeukoCt/Data_In/RAWVids_In` directory, performs
  flat line correction to remove illumination differences across each
  video frame, and equalize the tissue intensity between the video
  frames. The output is stored in `LeukoCt/Data_In/Reg`

###### INPUT:

- Videos in .tif or h5 .format, labelled in the scheme
  `[DATE]_P[PatientID]_V[VesselID]_[VideoID]_S[SliceStart-SliceEnd]_[a-z]`.
  Slice (= video frame) information and \[a-z\] (allows multiple vessels
  in one video) are optional.
- Dimensions of the videos, in the Excel file
  `LeukoCt/Data_In/OBMVideoProperties_In.xlsx`. Columns B-D specify the
  pixel/um value of each video pre-registration, at-registration, and
  post-registration respectively.

###### OPERATION:<br>

Update the subfolder names marked “\*\*\*\*\*” in the section
“User-Specified Variables”. Run script.

###### OUTPUT:<br>

Stored in `LeukoCt/Data_Out/Reg/Batch_FastStackReg_[TimeStamp]`.
Registered, flat-field corrected videos in at-registration and
post-registration dimensions, in .tif and .h5 format.  
<br> <br> \#### `STEP01_VesselROI.m` The script defines a ROI in the
input video with user assistance. It crops the video to exclude
far-from-ROI regions, and defines the time units for the skeleton
coordinate system.

This is the only script in the pipeline that requires user interaction
during its run.

###### INPUT:

- A registered video with .tif or .h5 extension, labelled in the Sample
  ID format. Registered videos for DEMO are stored in
  `LeukoCt/Data_In/REGVids_In`.
- Dimensions of the video, in the Excel file
  `LeukoCt/Data_In/OBMVideoProperties_In.xlsx`. Column E specifies the
  pixel/um value of the video in STEP01-07. Column F specifies the frame
  rate, in frames-per-second.
- Reference WBCC value in K/uL for each patient (Column B), identified
  by PatientID (Column A) from phlebotomy, in the Excel file
  `LeukoCt/Data_In/CBC_In.xlsx`.

###### OPERATION:

1.  Update the subfolder names marked “\*\*\*\*\*” in the section
    “User-Specified Variables”. Run script.
2.  User is asked to identify the vessel of interest in the video. Using
    the mouse, click anywhere on the vessel of interest, which starts a
    blue line, then end blue line draw anywhere with double click.
3.  After the software creates a vessel luminal mask, the user is asked
    to apply the “line eraser”. The line eraser (green line) breaks the
    mask where the drawn line crosses, and the software keeps the
    largest piece of mask after the break. Using the mouse, the user can
    apply the line eraser to remove undesirable areas, such as
    out-of-focus regions and where segments of the vessels intersect, or
    correct errors of the mask. The line eraser draw also terminates
    with a double click.  
4.  The software shows the updated vessel luminal mask overlaid with its
    skeleton, which colocalises with the flow axis except for the ends.
    User input is then requested to finalize the ROI. With the mouse,
    box in (in cyan) where the skeleton aligns well with the flow axis,
    starting with the end closer to the origin of the flow.
5.  The software completes the rest of the script.

###### OUTPUT:<br>

The script opens a folder of name `[SampleID]_PIPELINE_[TimeStamp]` in
`LeukoCt/Data_Out/Pipeline_MANUAL`. Output of this and subsequent steps
02-06 will be stored in its own subfolder, with a Workspace\_\*.mat file
containing the key outputs.

The key output variables of each step can be identified by variable
names starting with “o” (for output). These “o-variables” are highly
detailed; a more user-friendly, summarized output will be given in the
last step (`WBCConcStabilityAnalysis.m`; see below).

For STEP01, the key output is in the variable oROI. The saved image
`Img_PolygonROI*.tif` shows the vessel luminal mask, the skeleton (flow
axis) and the final ROI box drawn by the user. <br> <br> \####
`STEP02_SCS_VesselSz.m` \#### `STEP03_Velocity.m` \####
`STEP04_Velocity_Volume.m` \#### `STEP05_WBCCount.m` \####
`STEP06_WBCConc.m` \#### `STEP07_Report.m`

#### `STEP02_to_07_Batch_Pipeline.m`

- The script for STEP02 creates the skeleton coordinate system for the
  ROI and measure the vessel diameter.
- The script for STEP03 processes the space-time diagrams of each SCS
  unit, performs statistics on the slope angles and blood flow velocity.
- The script for STEP04 plots the blood flow velocity and calculates the
  blood flow volume rate and volume. It also flags SCS units with large
  uncertainties in these measurements.  
- The script for STEP05 counts leukocytes and performs all related
  tasks, including cell count window optimization.
- The script for STEP06 takes the results from STEP04 (blood flow
  volume) and STEP05 (leukocyte count) and estimates the leukocyte
  concentration.
- The script for STEP07 extracts the key outcomes from STEP01-07 and
  organizes them in an Excel Table format.

###### INPUT:<br>

(See OPERATION)

###### OPERATION:<br>

STEP02 to STEP07 can be run manually on individual videos, or
automatically and in batch mode, on multiple videos with their STEP01
`Workspace_STEP01Var_*.mat` output placed in the subfolder
`LeukoCt/Data_In/VesselROI_In`.

If run manually, no further input is necessary. Run the scripts in order
of the steps. Output is found in the same folder as STEP01 in
`LeukoCt/Data_Out/Pipeline_MANUAL`.

If run automatically and in batch mode, copy the STEP01
Workspace_STEP01Var\_\*.mat of all videos into the subfolder
`LeukoCt/Data_In/VesselROI_In`. Run the script
<u>STEP02_to_07_Batch_Pipeline.m</u>, which runs STEP02 to STEP07.
Output is found in `LeukoCt/Data_Out/Pipeline_BATCHAUTO`; each video’s
output is in the subfolder `[SampleID]_AUTO_PIPELINE_[TimeStamp]`.

###### OUTPUT:<br>

Key output variables:

- `STEP02_SCS_VesselSz.m`:
  - `oSCS` (skeleton coordinate system parameters)
  - `oRadius` (Vessel Lumen (L) Radius)
  - `oFocusQ` (ROI Contrast Map Output)
- `STEP03_Velocity.m`:
  - `oOrient` (Space-time diagram line slopes)
  - `oVelocity` (Flow Velocity)
- `STEP04_Velocity_Volume.m`:
  - `oVolume` (Flow Volume Rate (oVolume.VR) and Flow Volume
    (oVolume.V))
  - `oFlag_VctyEr` (SCS flagged for uncertain flow velocity and volume
    measurements)
- `STEP05_WBCCount.m`:
  - \`oCellCt (Time traces and cell count window parameters
    (oCellCt.Ct_mCell{1,4}))
- `STEP06_WBCConc.m`:
  - `oCellConc_mCell`

Other noteworthy outputs:

- `STEP02_SCS_VesselSz.m`:
  - `SkeletonCoordinateSystem_RefSlice*.tif` (SCS overlaid on ROI)
  - `Img_FocusDiameter_Illus_[TimeStamp].tif` (ROI Contrast Map)
  - `tMidPt_Vessel_Lumen_Diameter_um_Plot_*.tif` (Vessel luminal
    diameter along ROI)
- `STEP04_Velocity_Volume.m`:
  - `Plot_tMidPt_Velocity_um_s*.tif` (Flow velocity, in um/s)
  - `Plot_tMidPt_Volume_nL_min*.tif` (Flow volume rate, in nL/min)
- `STEP05_WBCCount.m`:
  - `PxWin_ImgMaskOverlay_WBCtoGateXYWinAreaRatio*.tif` (Final candidate
    cell count windows, of optimized size but different positions)
  - `PxWin_GMM_Plot_IPCellCt_WBCtoGateXYWinAreaRatio*.tif` (Time trace
    of leukocyte peaks for the candidate windows; “Optimized = 1” is the
    chosen window for leukocyte counting. Red dotted lines indicate the
    leukocyte-positive video frames from manual cell counting)
- `STEP07_Report.m`:
  - `[SampleID]_ReportDataLabel.csv` (important data from STEPS01-06 in
    spreadsheet format)
- `STEP02_to_07_Batch_Pipeline.m`:
  - `BatchRun_Data_w_Label*.csv` (same as `*ReportDataLabel.csv` from
    STEP07, but for multiple videos)
  - `BatchRun_Summary*.csv` (this spreadsheets indicates which video in
    the batch run fails to run properly, and at which step) <br> <br>

#### `WBCConcStabilityAnalysis.m`

This script combines data from the videos of the same ROI, and performs
stability analysis that determines if the leukocyte concentration has
stabilized sufficiently to be a robust predictor of the clinical WBCC
value.

###### INPUT:<br>

Pipeline outputs from STEPS00-07, with folder names
`[SampleID]_PIPELINE_[TimeStamp]` or
`[SampleID]_AUTO_PIPELINE_[TimeStamp]`. `LeukoCt/Data_In/Stability_In`
contain such outputs for DEMO.

###### OPERATION:<br>

Modify the subfolder names marked “\*\*\*\*\*” in the section
“User-Specified Variables”. Run script.

###### OUTPUT:<br>

Folder in `LeukoCt/Data_Out/Stability`: \*
`GroupID_AllVidData_WBCConcStability_Tbl_[TimeStamp].csv` summarizes the
stability analysis outcome for each ROI, coded as GroupID (format:
`[PPVV[a-z]]`, `PP` = PatientID, `VV` = VesselID, `[a-z]` is the same as
in SampleID). Column B is the vessel luminal diameter (um), Column C is
the stabilization time and Column D, the stabilization volume. -9999 for
the latter two values indicate the GroupID failed to achieve
concentration stability. \* Subfolder `Plot_FlowVolBin_CumWBCConc`
contains plots of cumulative leukocyte concentration for each GroupID.
\* `oASC` is the key output variable.
oASC.GroupID_AllVidData_gCell{\[GroupIDIdx\],2} contains the important
data from STEP01-07 for each ROI. The description of each column can be
found in the custom function
`LeukoCt/Dependencies/WBCConcStabilityAnalysis/fCombineVesselVidData`.

## Citation

Weighted statistics for slope and flow velocity calculations is from the
following reference:

Bevington, P. R., Data Reduction and Error Analysis for the Physical
Sciences, 336 pp., McGraw-Hill, 1969.

The following packages from [MATLAB File
Exchange](https://www.mathworks.com/matlabcentral/fileexchange/) were
used by the software:

1.  Manuel Guizar (2023). [Efficient subpixel image registration by
    cross-correlation](https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation),
    MATLAB Central File Exchange. Retrieved July 10, 2023.
2.  Dmitry Kaplan (2023). [Knee
    Point](https://www.mathworks.com/matlabcentral/fileexchange/35094-knee-point),
    MATLAB Central File Exchange. Retrieved June 25, 2023.
3.  Dirk-Jan Kroon (2022). [2D Line Curvature and
    Normals](https://www.mathworks.com/matlabcentral/fileexchange/32696-2d-line-curvature-and-normals),
    MATLAB Central File Exchange. Retrieved February 28, 2022.
4.  Oliver Woodford (2017). [real2rgb &
    colormaps](https://www.mathworks.com/matlabcentral/fileexchange/23342-real2rgb-colormaps),
    MATLAB Central File Exchange. Retrieved February 7, 2017.
5.  Zhaoxu Liu / slandarer (2024). [200
    colormap](https://www.mathworks.com/matlabcentral/fileexchange/120088-200-colormap),
    MATLAB Central File Exchange. Retrieved February 7, 2024.
6.  Sven Haase (2024). [weighted
    median](https://www.mathworks.com/matlabcentral/fileexchange/23077-weighted-median),
    MATLAB Central File Exchange. Retrieved February 23, 2024.
7.  Yair Altman (2024).
    [export_fig](https://github.com/altmany/export_fig/releases/tag/v3.40),
    GitHub. Retrieved June 17, 2024.
8.  Jos van der Geest (2024).
    [CATSTRUCT](https://www.mathworks.com/matlabcentral/fileexchange/7842-catstruct),
    MATLAB Central File Exchange. Retrieved February 6, 2024.
9.  Povilas Karvelis (2024). [daviolinplot - beautiful violin and
    raincloud
    plots](https://github.com/povilaskarvelis/DataViz/releases/tag/v3.2.4),
    GitHub. Retrieved March 21, 2024.
10. Chad Greene (2024).
    [Label](https://www.mathworks.com/matlabcentral/fileexchange/47421-label),
    MATLAB Central File Exchange. Retrieved January 14, 2024.
11. Adam Danz (2024).
    [labelpoints](https://www.mathworks.com/matlabcentral/fileexchange/46891-labelpoints),
    MATLAB Central File Exchange. Retrieved March 29, 2024.
