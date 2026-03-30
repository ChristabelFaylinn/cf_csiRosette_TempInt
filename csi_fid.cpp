//-----------------------------------------------------------------------------
// <copyright file="csi_fid_SAFETY.cpp" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
// <description>CSI FID sequence, no in-plane VoI selection</description>
//-----------------------------------------------------------------------------

#include "MrSpecAcq/cf_csiRosette_TempInt/csi_fid.h"

#include "MrImagingFW/libSeqUTIF/libsequt.h"
#include "MrImagingFW/libSeqUtilFW/SeqTrace.h"
#include "MrMeasSrv/MeasNuclei/IF/MeasKnownNuclei.h"
#include "MrMeasSrv/SeqFW/libGSL/libGSL.h"
#include "MrMeasSrv/SeqFW/libSSL/libSSL.h"
#include "MrMeasSrv/SeqIF/Sequence/Sequence.h"
#include "MrMeasSrv/SeqIF/Sequence/sequmsg.h"
#include "MrMeasSrv/SeqIF/csequence.h"
#include "MrMeasSrv/SeqIF/libRT/libRTDefines.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Filter/MrFilter.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/KSpace/MrKSpace.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/Physiology/MrPhysiology.h"
#include "MrSpecAcq/SpecUtils/EGA/SpecUtils.h"
#include "MrSpecAcq/cf_csiRosette_TempInt/csi_fid_UI.h"

// ------------------------------------------------------------------------------ 
// WIP MEM BLOCK includes                                                           
// ------------------------------------------------------------------------------ 
#include "MrProtSrv/Domain/CoreNative/MrWipMemBlock.h"

#include <algorithm>
#include <cmath>
#include <tuple>

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif

#ifdef SEQUENCE_CLASS_CSI_FID
SEQIF_DEFINE(SEQ_NAMESPACE::Csi_fid)
#endif

using namespace SEQ_NAMESPACE;

// Defining Global Variables
static sROT_MATRIX sNewRotationMatrix, sExistingMatrix;
static long /*lI,*/ lJ, lK, lL;
static int          n_ti = 2; // CF change to user input

static double ZRotMatrix[3][3];
static double dRotAngle;


using namespace std;

#pragma warning (push)
#pragma warning(disable : 4355)

#pragma warning(pop)

Csi_fid::~Csi_fid()
{
#ifdef WIN32
    if (m_pUI)
    {
        delete m_pUI;
        m_pUI = nullptr;
    }
#endif
}

NLSStatus Csi_fid::initialize(SeqLim& rSeqLim)
{
    if (SysProperties::isFreeDot())
    {
        return MRI_SEQ_SEQU_SEQU_SYSTEM_INCOMPATIBLE;
    }

    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    SpecSeq::initialize(rSeqLim);

    char adcName[15];
        for (int m=0;m<9;m++){
            sprintf(adcName,"ADC_%d",m);
            m_adc1[m].setIdent(adcName);
            }

    rSeqLim.setMyOrigFilename(__FILE__);
    rSeqLim.setSequenceOwner(SEQ_OWNER_SIEMENS);
    rSeqLim.setSequenceHintText("\n\
                Application: Spectroscopy \n\
                Basics: CSI, Spin-Echo \n\
                Build: " __DATE__ "   " __TIME__ "\n");

    rSeqLim.isSVSSequence(true);
    rSeqLim.isCSISequence(true);
    rSeqLim.setSequenceCard(SEQ::SEQUENCE_CARD_SPECTROSCOPY);
    rSeqLim.setSequenceType(SEQ::SEQUENCE_TYPE_CSI);
    rSeqLim.getSequenceType().setDisplayMode(SEQ::DM_OFF);
    SeqUT.setSequenceType(UT_SEQUENCE_TYPE_CSI_FID);

    // the system requirements: gradient power
    rSeqLim.setRequiredGradAmpl(15.0);      //  mT/m but we have 20 on TURBO gradients
    rSeqLim.setRequiredGradSlewRate(25.0);  // (mT/m)/ms, TURBO gradients
    rSeqLim.setGradients(SEQ::GRAD_NORMAL); //  SEQ::GRAD_FAST, SEQ::GRAD_WHISPER);

    SpecUtils::setForcePositioningOnNDIS(rSeqLim);

    // VectorSize of time domain signal
    // this is a spectroscopy specific variable
    rSeqLim.setVectorSize(128, 2048, SEQ::BASE2, 512);
    rSeqLim.setReadoutOSFactor(2.); // default
    rSeqLim.setRemoveOversampling(SEQ::YES, SEQ::NO);

    // Base matrix size of the image
    rSeqLim.setBaseResolution(8, 128, SEQ::INC_NORMAL, 48);
    rSeqLim.setPELines(8, 128, 1, 48);
    rSeqLim.setMaxPhaseResolution(2);

    rSeqLim.setfinalMatrixSizeRead(8, 128, SEQ::BASE2, 64);
    rSeqLim.setfinalMatrixSizePhase(8, 128, SEQ::BASE2, 64);

    // slices and partitions
    rSeqLim.setSlices(1, 1, 1, 1);


    rSeqLim.setContrasts(2, 2, 1, 2);

    // 3D CSI
    rSeqLim.setDimension(SEQ::DIM_2, SEQ::DIM_3); // order sensitive!, to generate a consistent 2D protocol, it is
                                                  // mandatory to set some default values in fSeqPrep()
    rSeqLim.setPartition(8, 16, SEQ::INC_NORMAL, 8);
    rSeqLim.setfinalMatrixSizeSlice(8, 16, SEQ::BASE2, 8);
    rSeqLim.setImagesPerSlab(8, 32, SEQ::BASE2, 8); // should be redundant!
    rSeqLim.setSlabThickness(200, 400);
    rSeqLim.set3DPartThickness(10, 50, 1, 50);
    rSeqLim.setMinSliceResolution(0.5); // why?

    // Bandwidth used for data acquisition
    rSeqLim.setBandWidth(0, 4200, 4200, 100, 4200); // CF bw of acquisition ??

    // Echo Time
    rSeqLim.setTE(0, 1400, 300000, 10, 15000);

    // Repetition Time
    rSeqLim.setTR(0, 200000, 30000000, 10000, 1500000);

    rSeqLim.setRepetitions(0, 7, 1, 0); // 2D max 8, 3D max 4 realized with GetLimitHandler
    rSeqLim.setRepetitionsDelayTime(0, 10000000, 100000, 0);

    // CSI FOV
    rSeqLim.setReadoutFOV(200, 400, 10, 240);
    rSeqLim.setPhaseFOV(200, 400, 10, 240);
    rSeqLim.setSliceThickness(5.0, 80.0, 1, 25.0);
    rSeqLim.getSliceDistanceFactor().setDisplayMode(SEQ::DM_OFF);

    // VOI definition
    // these are spectroscopy specific variables
    rSeqLim.setVoIPosCor(-150, 150, .1, 0);
    rSeqLim.setVoIPosSag(-150, 150, .1, 0);
    rSeqLim.setVoIPosTra(-150, 150, .1, 0);

    rSeqLim.setVoISizePhase(200, 400, 10, 400);
    rSeqLim.setVoISizeReadout(400, 400, 10, 400);
    rSeqLim.setVoISizeSlice(80, 200, 10, 80);

    // RF
    rSeqLim.setFlipAngle(0.000, 180.000, 1.000, 90.000);
    rSeqLim.setExtSrfFilename("%MEASCONST%/extrf/extrf_spec.pls");

    // phase encoding type
    // this is a spectroscopy specific variable
    rSeqLim.setPhaseEncodingType(
        SEQ::PHASE_ENCODING_FULL, SEQ::PHASE_ENCODING_WEIGHTED, SEQ::PHASE_ENCODING_ELLIPTICAL);

    // water suppression (default setting )
    rSeqLim.setFatWaterContrast(
        FatWaterContrast_WaterSaturation,
        FatWaterContrast_Standard,
        FatWaterContrast_WaterSuppressionWeak,
        FatWaterContrast_WaterSuppresionRFOff);

    // Preparation Pulses
    // this is a spectroscopy specific variable
    rSeqLim.setPreparingScans(0, 16, 1, 4);

    // Acquisition delay
    // this is a spectroscopy specific variable
    rSeqLim.setAcquisitionDelay(1400, 4500, 100, 2300);

    // Averages / Repetitions
    rSeqLim.setAverages(1, 4096, 1, 2);

    // Bandwidth of Water Excitation Pulses
    rSeqLim.setRfBandwidth(20, 80, 5, 35);

    // Loop control
    rSeqLim.setAveragingMode(SEQ::INNER_LOOP, SEQ::OUTER_LOOP);

    // delta frequency range of TX/RX
    rSeqLim.setDeltaFrequency(0, -10000, 10000, 1, 0);

    // delta frequency range of TX nucleus
    rSeqLim.setDeltaFrequency(1, -10000, 10000, 1, 0);

    // spectroscopy sequences do NOT use the raw data filters of the imaging sequences
    rSeqLim.setFilterType(SEQ::HAMMING, SEQ::PRESCAN_NORMALIZE);

    // NOE settings
    rSeqLim.setNOEs(1, 10, 1, 1); // up to 10 NOE pulses
    rSeqLim.setNOEType(SEQ::NOE_NONE, SEQ::NOE_RECTANGULAR);
    rSeqLim.setNOEDuration(1000, 90000, 1000, 5000);
    rSeqLim.setNOEDelay(1000, 100000, 1000, 10000);
    rSeqLim.setNOEFlipAngle(0, 180, 1, 90);

    // DC settings
    rSeqLim.setDecouplingType(SEQ::DECOUPLING_NONE, SEQ::DECOUPLING_WALTZ4, SEQ::DECOUPLING_CW);
    rSeqLim.setDecouplingFlipAngle(0, 360, 1, 180);
    rSeqLim.setDecouplingDuration(500, 2000, 500, 1000); // the duration of the 180 deg. WALTZ pulse
    rSeqLim.setDecouplingTotalDuration(10, 100, 10, 50); // in % of acquisition time
    rSeqLim.setDecouplingPause(0, 200, 10, 100);         // percentage value

    // Physiologic measurements
    rSeqLim.addPhysioMode(SEQ::SIGNAL_NONE, SEQ::METHOD_NONE);
    rSeqLim.addPhysioMode(SEQ::SIGNAL_CARDIAC, SEQ::METHOD_TRIGGERING);
    rSeqLim.addPhysioMode(SEQ::SIGNAL_RESPIRATION, SEQ::METHOD_TRIGGERING);
    rSeqLim.setPhases(1, 1, 1, 1); //  phases not supported in loop structure

    // default adjust procedures
    rSeqLim.setAdjWatSup(SEQ::ENABLE, SEQ::DISABLE);

    // new with VA21A
    rSeqLim.setSliceSelectDeltaFrequency(-10.0, 10.0, 0.1, 0.0); // suitable range for X-nuclei

    // the coil combine mode, we don't want to have; adaptive coil combine is not possible, because then save uncombined
    // is not allowed
    ParLimOption<SEQ::CoilCombineMode>& coilCombine = rSeqLim.getCoilCombineMode();
    coilCombine.set(SEQ::COILCOMBINE_SUM_OF_SQUARES);
    coilCombine.setDisplayMode(SEQ::DM_OFF);

    // Data Receive & Image calculation
    rSeqLim.setICEProgramFilename("%SiemensIceProgs%\\IcePrgSpectroscopy");

    // Deactivate channel reduction (software matrix mode)
    rSeqLim.getChannelDiscardMode().setDisplayMode(SEQ::DM_OFF);
    rSeqLim.getChannelMixingMode().setDisplayMode(SEQ::DM_OFF);

    // support of the other nuclei
    rSeqLim.setSupportedNuclei(NUCLEI_ALL.get().c_str());

#ifdef WIN32
    //  ----------------------------------------------------------------------
    //  Instantiate of UI class
    //  ----------------------------------------------------------------------
    createUI();

    //  ----------------------------------------------------------------------
    //  Declaration of pointer to UI parameter classes
    //  ----------------------------------------------------------------------
    lStatus = m_pUI->registerUI(rSeqLim);

    if (NLS_SEVERITY(lStatus) != NLS_SUCCESS)
    {
        SEQ_TRACE_ERROR.print("Initialization of UI failed : 0x%lx", lStatus);
        return lStatus;
    }
#endif
    return lStatus;
}

NLSStatus Csi_fid::prepare(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    MrProtocolData::SeqExpoRFInfo rfInfoMainNucPerScan;
    MrProtocolData::SeqExpoRFInfo rfInfoTxNucPerScan;

    // long sp1dur, sp2dur;
    // long encodTotaldur; //ss
    double /*sp1ampl, sp2ampl, */ dmin, dmax;

    if (rMrProt.preScanNormalizeFilter().getucOn())
        rMrProt.preScanNormalizeFilter().setucStoreCXIma(true); // exception, change a protocol during prepare() !!!

    if (rMrProt.kSpace().dimension() == SEQ::DIM_2)
    {
        rMrProt.getsSpecPara().setlFinalMatrixSizeSlice(1); // finalMatrixSizeSlice is 1 in 2D case
        rMrProt.kSpace().setlPartitions(1);                 // partitions is 1 in 2D case
        rMrProt.getsSpecPara().getsVoI().setdThickness(
            rMrProt.sliceSeries()[0].thickness()); // VoiSizeSlice = SliceThickness
    }

    // Get the current nucleus from the protocol
    MeasNucleus  mainNucleus(rMrProt.getsTXSPEC().getasNucleusInfo()[0].gettNucleus().c_str());
    const double dLarmorConst = mainNucleus.getLarmorConst();
    MeasNucleus  txNucleus(rMrProt.getsTXSPEC().getasNucleusInfo()[1].gettNucleus().c_str());

    // do some basic checking first
    if (txNucleus == NUCLEUS_NONE)
    {
        if (rMrProt.getsSpecPara().getlNOEType() != SEQ::NOE_NONE)
        {
            SEQ_TRACE_DEBUG.print("NOE requires the Tx nucleus to be defined.");
            return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;
        }

        if (rMrProt.getsSpecPara().getlDecouplingType() != SEQ::DECOUPLING_NONE)
        {
            SEQ_TRACE_DEBUG.print("Decoupling requires the Tx nucleus to be defined.");
            return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;
        }
    }
    else // a TX nucleus is selected
    {
        if (txNucleus != NUCLEUS_1H)
        {
            SEQ_TRACE_DEBUG.print("The Tx nucleus must be 1H");
            return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;
        }
    }

    // the coil combine mode, we don't want to have; adaptive coil combine is not possible, because then save uncombined
    // is not allowed
    if (rMrProt.coilCombineMode() != SEQ::COILCOMBINE_SUM_OF_SQUARES)
        return MRI_SEQ_SEQU_ERROR;

    // compute VoI and FoV orientation and position
    if (!(m_voi.prep(rMrProt, rSeqLim, rMrProt.getsSpecPara().getsVoI(), 0)))
        return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;

    // in the case of multi-slice CSI we would need to prepare more than the one slice
    // in this case, the function fSUPrepSlicePosArray() might be helpful
    if (!(m_fov.prep(rMrProt, rSeqLim, rMrProt.sliceSeries()[0], 0)))
        return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;

    // compute general limitations
    const long   lRampTime_us = 200;                                             // us
    const double dGradMaxAmpl = SysProperties::getGradMaxAmpl(SEQ::GRAD_NORMAL); // max gradient amplitude

    // spoiler sp1:  sp1dur * sp1ampl = 14000 us*mT/m
    // sp1dur = 700;	
    // sp1ampl = 14000./sp1dur;

	// spoiler sp2:  sp2dur * sp2ampl = 14000 us*mT/m... but does this need to be 46000us*mT/m? /lb
    // sp2dur = 700;	
    // sp2ampl = 14000./sp2dur;

    // set product sequence default values
    const long   l_ExcDuration = 1280;
    // const long   l_RefDuration = 5200; //lb this duration is specific to MAO pulse, hard-coded value
    const long   encoddur      = 500;
    // const long   encodTotaldur = encoddur+lRampTime_us;
    const double refoc         = 0.515;
    const double dVoxelShiftR  = 0.5;
    const double dVoxelShiftP  = 0.5;
    const double dVoxelShiftS  = 0.5;

    // Prepare Osc. Bit
    m_osc1.setCode(SYNCCODE_OSC0);
    m_osc1.setIdent(RTEIDENT_Osc0);
    m_osc1.setDuration(10); // 10 see ...\MrMeasSrv/SeqIF/csequence.h

    //////////////////////////////////////////////
    // excitation pulse (SLICE direction)
    ///////////////////////////////////// ////////

    // excitation and refocussing pulse (SLICE direction)
    m_rf_exc.setTypeExcitation();
    m_rf_exc.setDuration(l_ExcDuration);
    m_rf_exc.setFlipAngle(rMrProt.flipAngle());
    m_rf_exc.setInitialPhase(0);
    m_rf_exc.setFamilyName("hsinc_512_3550");
    m_rf_exc.setNucleus(mainNucleus);
    m_rf_exc.setThickness(rMrProt.getsSpecPara().getsVoI().getdThickness());

    if (!(m_rf_exc.prepExternal(rMrProt, rSeqExpo)))
        return m_rf_exc.getNLSStatus();

    rfInfoMainNucPerScan += m_rf_exc.getRFInfo();

    // m_rf_ref.setTypeRefocussing();
    // m_rf_ref.setDuration(l_RefDuration);
    // m_rf_ref.setFlipAngle(180);
    // m_rf_ref.setInitialPhase(0);
    // m_rf_ref.setFamilyName("mao_400_4");
    // m_rf_ref.setNucleus(mainNucleus);
    // m_rf_ref.setThickness(rMrProt.getsSpecPara().getsVoI().getdThickness());

    // if (!(m_rf_ref.prepExternal(rMrProt, rSeqExpo)))
    //     return m_rf_ref.getNLSStatus();

    // rfInfoMainNucPerScan += m_rf_ref.getRFInfo();
    // computation of the frequency offset which defines the voxel position
    double dFrequency_Hz        = 0.0;
    double dFrequencyWithCSD_Hz = 0.0;

    std::tie(dFrequency_Hz, dFrequencyWithCSD_Hz) = SpecUtils::calcFrequency(
        m_rf_exc.getGSAmplitude(),
        dLarmorConst,
        m_fov.getSliceShift(),
        rMrProt.getsTXSPEC().getasNucleusInfo()[0].getlFrequency(),
        rMrProt.getsSpecPara().getdDeltaFrequency());

    SpecUtils::setExcFreqAndPhase(
        m_ph_s_exc, m_ph_n_exc, dFrequencyWithCSD_Hz, m_rf_exc.getDuration(), m_rf_exc.getAsymmetry());

    // double dFrequency_ref_Hz        = 0.0;
    // double dFrequencyWithCSD_ref_Hz = 0.0;

    // std::tie(dFrequency_ref_Hz, dFrequencyWithCSD_ref_Hz) = SpecUtils::calcFrequency(
    //     m_rf_ref.getGSAmplitude(),
    //     dLarmorConst,
    //     m_fov.getSliceShift(),
    //     rMrProt.getsTXSPEC().getasNucleusInfo()[0].getlFrequency(),
    //     rMrProt.getsSpecPara().getdDeltaFrequency());

    // SpecUtils::setExcFreqAndPhase(
    //     m_ph_s_ref, m_ph_n_ref, dFrequencyWithCSD_ref_Hz, m_rf_ref.getDuration(), m_rf_ref.getAsymmetry());


    // prepare the NOE pulse(s)
    long noe_dur = 0;
    if (rMrProt.getsSpecPara().getlNOEType() == SEQ::NOE_RECTANGULAR)
    {
        m_rf_rect_noe.setNucleus(txNucleus);
        long lDuration = rMrProt.getsSpecPara().getlNOEDuration();
        lDuration = (lDuration % 100) ? (lDuration / 100) * 100 : lDuration; // rounding down to a multiple of 100 us
        m_rf_rect_noe.setDuration(lDuration);
        m_rf_rect_noe.setSamples(lDuration / 100);
        m_rf_rect_noe.setFlipAngle(rMrProt.getsSpecPara().getdNOEFlipAngle());
        m_rf_rect_noe.setInitialPhase(0.0);  // not used
        m_rf_rect_noe.setThickness(10000.0); // not used
        if (!m_rf_rect_noe.prepRect(rMrProt, rSeqExpo))
            return m_rf_rect_noe.getNLSStatus();

        noe_dur = rMrProt.getsSpecPara().getlNOEs() * (lDuration + rMrProt.getsSpecPara().getlNOEDelay());
        rfInfoTxNucPerScan += m_rf_rect_noe.getRFInfo() * rMrProt.getsSpecPara().getlNOEs();

        m_ph_s_noe.setNucleus(txNucleus);
        m_ph_n_noe.setNucleus(txNucleus);

        m_ph_s_noe.setFrequency((long)(rMrProt.getsTXSPEC().getasNucleusInfo()[1].getlDeltaFrequency()));
        m_ph_n_noe.setFrequency(0L);
        m_ph_s_noe.setPhase(0);
        m_ph_n_noe.setPhase(0);
    }

    //  Prepare the readout event
    // m_adc1.setColumns(rMrProt.getsSpecPara().getlVectorSize()); // takes OS into account
    // m_adc1.setDwellTime(
        // 10 * std::lround(1E8 / rMrProt.bandWidth(rSeqLim.getReadoutOSFactor())[0])); // in ns, takes OS into account
    m_lSamplesBeforeEcho = 0; // we are always behind the echo in a FID sequence

    // Additions for compatibility with the functor-based Ice program
    // m_adc1.getMDH().setKSpaceCentreColumn(0);
    // m_adc1.getMDH().setKSpaceCentreLineNo(
    //     (uint16_t)rMrProt.kSpace().phaseEncodingLines() / 2); // necessary to reduce warnings from Ice program
    // m_adc1.getMDH().setKSpaceCentrePartitionNo(
    //     (uint16_t)rMrProt.kSpace().getlPartitions() / 2); // necessary to reduce warnings from Ice program

    // In place of IceProgramPara - do not include oversampling
    // m_adc1.getMDH().setPreCutOff(0);
    // m_adc1.getMDH().setPostCutOff(0);

    int Nx = rMrProt.kSpace().baseResolution();
    int ADCdwell=10000;
    // Initialize some values:
    double sw = rMrProt.bandWidth(rSeqLim.getReadoutOSFactor())[0]
                / n_ti; // CF This is hardcoded in the "initialize" section to be 2500 Hz.
    double fov = rMrProt.sliceSeries()[0].readoutFOV(); // Get the FOV
    double Y   = 42.577E6;                              // Gamma in Hz/T
    // long n_ti = 1;  //  Number of temporal interleaves (Fixed to 1 and not used.  For future use);
    double delta_t         = 0.00001; // gradient raster time in seconds.  0.00001 sec = 10 microsec
    long   loops_per_block = 85;      // Number of loops around a rosette petal (i.e. FID points) per ADC/rosette block.
    long   spectral_points = 765;     // Total number of loops (i.e. FID points) in all ADC/rosette blocks together.
    long   grad_sample_points;        // Total number of gradient sample points in a single ADC/rosette block bl

    // Calculated parameters
    double w        = M_PI * sw;                         // Rosette frequency omega (w = w1 = w2).
    double Nsh      = ceil((M_PI * Nx) / 2);             // Number of rosette petals (shots) Nsh.
    int    n_blocks = spectral_points / loops_per_block; // Number of ADC blocks to stitch together.
    double kmax     = Nx * 1000.0 / (2.0 * fov);         // kmax in [m^-1]
    // double delta_t_spectral = 1/sw; // NOT USED.

    // Rounding values of grad_sample_points
    // grad_sample_points=loops_per_block*ceil(M_PI/(w*delta_t));
    // grad_sample_points= ceil(grad_sample_points);
    grad_sample_points = 3995; // number of grad sample points is hard-coded for now.  Uncomment two above lines to
                               // calculate dynamically.

    double points_per_loop
        = (1 / static_cast<double>(sw))
          / (static_cast<double>(ADCdwell)
             / 2e9); // Number of ADC points per rosette loop: // CF ISNT THIS SUPPOSED TO JUST BE /1e9???
    // double dt_adc_max = (1/(static_cast<double>(sw)*(1+(M_PI*static_cast<double>(Nx)/2))));  //The max value of
    // dt_adc.  NOT USED.

    // CF  delay_3 to ensure that it is a multiple of 10ns

    SEQ_TRACE_WARN.print("sw before rounding: %d", loops_per_block / (delta_t * grad_sample_points ));

    w     = (M_PI * loops_per_block / (delta_t * grad_sample_points));
    dummy_delay_3 = M_PI * 1000000 / (n_ti * w);

    if ((dummy_delay_3 % 10) != 0)
    {
        std::cout << "Our original dummy value for delay_3 is: " << dummy_delay_3 << std::endl;
       /* float w_f
            = (M_PI * loops_per_block
               / (delta_t * grad_sample_points));*/ // new w to accomodate rounded grad_sample_points values
        //std::cout << "Our dummy value for w is: " << w << std::endl;
        //std::cout << "Our dummy value for w_f is: " << w_f << std::endl;
        //float alpha_here = w_f / w; // our scaling factor which we need to keep, CF when rounding happened in between
        //std::cout << "Our dummy value for alpha is: " << alpha_here << std::endl;
        dummy_delay_3 = (10 - ((dummy_delay_3) % 10)) + dummy_delay_3; // Round our delay_3 to a multiple of 10
        //std::cout << "Our new dummy value for delay_3 is: " << dummy_delay_3 << std::endl;
        w = (M_PI * 1000000) / (n_ti * dummy_delay_3); // redefine w_f to accomodate rounded delay_3 value
        //w   = w_f / alpha_here;                          // Return w with the conversion factor
        std::cout << "Our rounded w is: " << w << std::endl;
        grad_sample_points = loops_per_block * ceil(M_PI / (w * delta_t));
        grad_sample_points = ceil(grad_sample_points);
        points_per_loop    = (grad_sample_points / spectral_points) * 10; // check if it divi
    }

   

    // Recalculate omega (w) and sw based on above
    w  = M_PI * loops_per_block / (delta_t * grad_sample_points);
    sw = w / M_PI;

    SEQ_TRACE_WARN.print("grad_sample_points after rounding: %d", grad_sample_points);
    SEQ_TRACE_WARN.print("sw after rounding: %d", sw);

 /*   int Nx = rMrProt.kSpace().baseResolution();
        int ADCdwell;*/
        for(int m=0;m<9;m++)
        	{   //if(Nx==64)
        //    { 
                ADCdwell=10000;
                m_adc1[m].prep (grad_sample_points, static_cast<int32_t>(ADCdwell)); // CF change 4026 to generalize with n_ti and sw
        //    }
        //else if(Nx==48)
        //      {  
        //		  ADCdwell=10000;
        //		  m_adc1[m].prep (4032, static_cast<int32_t>(ADCdwell));
        //      }
        //else 
        //	{  
        //		ADCdwell=16000;
        //		m_adc1[m].prep (3360, static_cast<int32_t>(ADCdwell)); }
        }

    m_ph_s_adc.setFrequency((long)(rMrProt.getsTXSPEC().getasNucleusInfo()[0].getlDeltaFrequency()));
    m_ph_n_adc.setFrequency(0L);
    m_ph_s_adc.setPhase(0);
    m_ph_n_adc.setPhase(0);

    // prepare the DC pulse(s), after preparing the readout event
    m_txrfon.setCode(SYNCCODE_TXRFON_DC);
    m_txrfon.setDuration(m_adc1[0].getRoundedDuration(GRAD_RASTER_TIME) + GRAD_RASTER_TIME * 4);

    // prepare the DC pulse(s), after preparing the readout event
    long lDecTotDurInPerCent = rMrProt.getsSpecPara().getlDecouplingTotalDuration();
    if (lDecTotDurInPerCent < 10)
        lDecTotDurInPerCent = 10; // minimum in seq lim = 10
    if (lDecTotDurInPerCent > 100)
        lDecTotDurInPerCent = 100;

    m_ln_dc_pulses = 0;
    if (rMrProt.getsSpecPara().getlDecouplingType() == SEQ::DECOUPLING_CW)
    {
        m_rf_dc.setNucleus(txNucleus);
        // make 49.5 ms DC pulses with 500 us interrupts
        if ((m_ln_dc_pulses = (long)(m_adc1[0].getRoundedDuration() * lDecTotDurInPerCent * 0.01 / 50000.)) > 0)
        {
            // make 49.5 ms DC pulses with 500 us interrupts
            rMrProt.utilityParameter()[1] = m_ln_dc_pulses * 50; // DC total duration in ms
            m_rf_dc.setDuration(49500);
            m_rf_dc.setSamples(495);
        }
        else
        {
            // DC duration shorter than 50 ms, so, prepare a shorter pulse
            m_ln_dc_pulses = 1;
            m_rf_dc.setDuration((((long)(m_adc1[0].getRoundedDuration() * lDecTotDurInPerCent * 0.01)) / 100) * 100);
            m_rf_dc.setSamples(m_rf_dc.getDuration() / 100);
            rMrProt.utilityParameter()[1] = (long)(m_rf_dc.getDuration() * .001); // us -> ms
        }
        m_rf_dc.setFlipAngle(
            m_rf_dc.getDuration() * .0001
            * rMrProt.getsSpecPara().getdDecouplingFlipAngle()); // the flip-angle per 10 ms
        m_rf_dc.setInitialPhase(0.0);                            // not used
        m_rf_dc.setThickness(10000.0);                           // not used
        if (!m_rf_dc.prepRect(rMrProt, rSeqExpo))
            return m_rf_dc.getNLSStatus();

        rfInfoTxNucPerScan += m_rf_dc.getRFInfo() * m_ln_dc_pulses;

        m_ph_s_dc.setNucleus(txNucleus);
        m_ph_n_dc.setNucleus(txNucleus);
        m_ph_s_dc.setFrequency((long)(rMrProt.getsTXSPEC().getasNucleusInfo()[1].getlDeltaFrequency()));
        m_ph_n_dc.setFrequency(0L);
        m_ph_s_dc.setPhase(0);
        m_ph_n_dc.setPhase(0);
    }
    else if (rMrProt.getsSpecPara().getlDecouplingType() == SEQ::DECOUPLING_WALTZ4) // WALTZ_4 DC
    {
        long i, ii, j;

        const long sample_dur = 50;             // us
        m_lPulsePauseRequired = sample_dur * 5; // must be multiple of sample_dur and >= 150 us

        const long MAX_PULSE_DUR = 50000; // us, MAX_PULSE_DUR / sample_dur < WALTZ_MAX_NO_SAMPLES
        long       n_samples_180 = rMrProt.getsSpecPara().getlDecouplingDuration() / sample_dur;
        long       n_samples_90  = n_samples_180 / 2;
        long       n_samples_270 = n_samples_180 + n_samples_90;
        long       n_samples_pause
            = (long)((double)n_samples_180 * .01 * (double)rMrProt.getsSpecPara().getlDecouplingPauseFract());

        // rSeqLim.setDecouplingDuration( 500, 2000, 500, 1000); // the duration of the 180 deg. WALTZ pulse
        // rSeqLim.setDecouplingPause( 0, 200, 10, 100 ); // percentage value
        // ->
        // min pulse dur = 500, max. pulse dur = 2000
        // min pause dur = 250, max. pause dur = 4000
        // -> 3 pulse dur.,  min. = 1750, max. = 18000
        // -> WALTZ cycle dur., min = 7000, max. = 72000 (the latter is the true maximum pulse duration)

        long n_waltz_cycles
            = (MAX_PULSE_DUR / sample_dur) / (4 * (n_samples_180 + n_samples_90 + n_samples_270 + 3 * n_samples_pause));
        if (n_waltz_cycles < 1) // this might be the case since MAX_PULSE_DUR < max. WALTZ cycle dur.
            n_waltz_cycles = 1;

        // compute the samples of the arbitrary pulse realizing n_waltz_cycles
        long   pulse_dur = 0;
        double ampl      = 0;
        for (i = 0; i < n_waltz_cycles; i++)
        {
            for (ii = 0; ii < 2; ii++)
            {
                for (j = 0; j < n_samples_90; j++)
                { // 90 deg.
                    m_samples[pulse_dur].flAbs = 1.0;
                    ampl += 1.0;
                    m_samples[pulse_dur].flPha = 0;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_pause; j++)
                { // pause
                    m_samples[pulse_dur].flAbs = 0;
                    m_samples[pulse_dur].flPha = 0;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_180; j++)
                { // 180 deg.
                    m_samples[pulse_dur].flAbs = 1.0;
                    ampl += 1.0;
                    m_samples[pulse_dur].flPha = (float)M_PI;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_pause; j++)
                { // pause
                    m_samples[pulse_dur].flAbs = 0;
                    m_samples[pulse_dur].flPha = 0;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_270; j++)
                { // 270 deg.
                    m_samples[pulse_dur].flAbs = 1.0;
                    ampl += 1.0;
                    m_samples[pulse_dur].flPha = 0;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_pause; j++)
                { // pause
                    m_samples[pulse_dur].flAbs = 0;
                    m_samples[pulse_dur].flPha = 0;
                    pulse_dur++;
                }
            }
            for (ii = 0; ii < 2; ii++)
            {
                for (j = 0; j < n_samples_90; j++)
                { // 90 deg.
                    m_samples[pulse_dur].flAbs = 1.0;
                    ampl += 1.0;
                    m_samples[pulse_dur].flPha = (float)M_PI;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_pause; j++)
                { // pause
                    m_samples[pulse_dur].flAbs = 0;
                    m_samples[pulse_dur].flPha = 0;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_180; j++)
                { // 180 deg.
                    m_samples[pulse_dur].flAbs = 1.0;
                    ampl += 1.0;
                    m_samples[pulse_dur].flPha = 0;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_pause; j++)
                { // pause
                    m_samples[pulse_dur].flAbs = 0;
                    m_samples[pulse_dur].flPha = 0;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_270; j++)
                { // 270 deg.
                    m_samples[pulse_dur].flAbs = 1.0;
                    ampl += 1.0;
                    m_samples[pulse_dur].flPha = (float)M_PI;
                    pulse_dur++;
                }
                for (j = 0; j < n_samples_pause; j++)
                { // pause
                    m_samples[pulse_dur].flAbs = 0;
                    m_samples[pulse_dur].flPha = 0;
                    pulse_dur++;
                }
            }
        }

        pulse_dur
            = pulse_dur
              - m_lPulsePauseRequired / sample_dur; // last WALTZPulsDelay is reduced of slPulsePause to compensate
                                                    // necessary PulsePause between  single arb. pulse shapes

        // prepare the pulse
        m_rf_waltz.setNucleus(txNucleus);
        m_rf_waltz.setSamples(pulse_dur);
        m_rf_waltz.setDuration(pulse_dur * sample_dur);
        m_rf_waltz.setFlipAngle(n_waltz_cycles * 4.0 * 3.0 * rMrProt.getsSpecPara().getdDecouplingFlipAngle());
        m_rf_waltz.setInitialPhase(0.0);  // not used
        m_rf_waltz.setThickness(10000.0); // not used
        if (!(m_rf_waltz.prepArbitrary(rMrProt, rSeqExpo, m_samples.data(), ampl)))
            return m_rf_waltz.getNLSStatus();

        if ((m_adc1[0].getRoundedDuration() / (pulse_dur * sample_dur + m_lPulsePauseRequired))
            < 1) // if the acquisition window is too short for even 1 WALTZ cycle
        {
            SEQ_TRACE_ERROR.print(
                "The acquisition window of %ld us duration does not allow a single WALTZ4 cycle of %ld us duration",
                m_adc1[0].getRoundedDuration(),
                pulse_dur * sample_dur + m_lPulsePauseRequired);
            return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;
        }

        // find out how many of these arbitrary pulses are needed to fill the acquisition window
        m_ln_dc_pulses
            = (long)(m_adc1[0].getRoundedDuration() * lDecTotDurInPerCent * .01 / (double)(pulse_dur * sample_dur + m_lPulsePauseRequired));

        if (m_ln_dc_pulses < 1)
            m_ln_dc_pulses = 1;

        // compute the total decoupling duration in ms
        rMrProt.utilityParameter()[1]
            = (long)(.001 * m_ln_dc_pulses * (pulse_dur * sample_dur + m_lPulsePauseRequired));

        // and prepare the trivial FREQ/PHASE event
        m_ph_s_dc.setNucleus(txNucleus);
        m_ph_n_dc.setNucleus(txNucleus);
        m_ph_s_dc.setFrequency((long)(rMrProt.getsTXSPEC().getasNucleusInfo()[1].getlDeltaFrequency()));
        m_ph_n_dc.setFrequency(0L);
        m_ph_s_dc.setPhase(0);
        m_ph_n_dc.setPhase(0);
        rfInfoTxNucPerScan += m_rf_waltz.getRFInfo() * m_ln_dc_pulses;
    }

    // compute grid position
    m_d_read_pos  = m_fov.getSliceOffCenterRO();
    m_d_phase_pos = m_fov.getSliceOffCenterPE();
    m_d_slice_pos = m_fov.getSliceShift();

    long n = 0;
    if (!((n = rMrProt.getsSpecPara().getlFinalMatrixSizeRead()) % 2)) // this is currently always the case
        m_d_read_pos += dVoxelShiftR * rMrProt.sliceSeries().front().readoutFOV() / static_cast<double>(n);

    if (!((n = rMrProt.getsSpecPara().getlFinalMatrixSizePhase()) % 2)) // this is currently always the case
        m_d_phase_pos += dVoxelShiftP * rMrProt.sliceSeries().front().phaseFOV() / static_cast<double>(n);

    if (!((n = rMrProt.getsSpecPara().getlFinalMatrixSizeSlice()) % 2)) // this is currently always the case
        m_d_slice_pos += dVoxelShiftS * rMrProt.sliceSeries().front().thickness() / static_cast<double>(n);

    // Prepare the gradient pulse structures

    // Prepare Rosette gradients
    // PARAMETERS FOR ROSETTE GRADIENTS

  //  //Initialize some values:
		//double sw = rMrProt.bandWidth(rSeqLim.getReadoutOSFactor())[0] / n_ti;  // CF This is hardcoded in the "initialize" section to be 2500 Hz.
		//double fov = rMrProt.sliceSeries()[0].readoutFOV();  //Get the FOV
		//double Y = 42.577E6; // Gamma in Hz/T
		//// long n_ti = 1;  //  Number of temporal interleaves (Fixed to 1 and not used.  For future use);
		//double delta_t = 0.00001; // gradient raster time in seconds.  0.00001 sec = 10 microsec 
		//long loops_per_block = 85;  //Number of loops around a rosette petal (i.e. FID points) per ADC/rosette block.
  //      long spectral_points = 765;  //Total number of loops (i.e. FID points) in all ADC/rosette blocks together.
		//long grad_sample_points;  //Total number of gradient sample points in a single ADC/rosette block bl
		//
		//// Calculated parameters
		//double w = M_PI*sw;  //Rosette frequency omega (w = w1 = w2).
		//double Nsh =  ceil((M_PI*Nx)/2);  //Number of rosette petals (shots) Nsh.
		//int n_blocks = spectral_points/loops_per_block;  //Number of ADC blocks to stitch together.
		//double kmax = Nx*1000.0/(2.0*fov); //kmax in [m^-1]
		//// double delta_t_spectral = 1/sw; // NOT USED.
  //      
		//// Rounding values of grad_sample_points 
  //      //grad_sample_points=loops_per_block*ceil(M_PI/(w*delta_t));
		////grad_sample_points= ceil(grad_sample_points);
  //      grad_sample_points = 3995;  //number of grad sample points is hard-coded for now.  Uncomment two above lines to calculate dynamically. 

		//double points_per_loop = (1/static_cast<double>(sw))/(static_cast<double>(ADCdwell)/2e9);  //Number of ADC points per rosette loop: // CF ISNT THIS SUPPOSED TO JUST BE /1e9???
		//// double dt_adc_max = (1/(static_cast<double>(sw)*(1+(M_PI*static_cast<double>(Nx)/2))));  //The max value of dt_adc.  NOT USED.

  //      // CF  delay_3 to ensure that it is a multiple of 10ns

  //      dummy_delay_3 = M_PI * 1000000 / (n_ti * w);
  //      
  //      if ((dummy_delay_3 % 10) != 0)
  //      {
  //          std::cout << "Our original dummy value for delay_3 is: " << dummy_delay_3 << std::endl;
  //          float w_f = (M_PI * loops_per_block / (delta_t * grad_sample_points)); // new w to accomodate rounded grad_sample_points values
  //          std::cout << "Our dummy value for w is: " << w << std::endl;
  //          std::cout << "Our dummy value for w_f is: " << w_f << std::endl;
  //          float alpha_here = w_f / w; // our scaling factor which we need to keep, CF when rounding happened in between 
  //          std::cout << "Our dummy value for alpha is: " << alpha_here << std::endl;
  //          dummy_delay_3 = (10 - ((dummy_delay_3) % 10)) + dummy_delay_3; // Round our delay_3 to a multiple of 10
  //          std::cout << "Our new dummy value for delay_3 is: " << dummy_delay_3 << std::endl;
  //          w_f = (M_PI * 1000000) / (n_ti * dummy_delay_3); // redefine w_f to accomodate rounded delay_3 value
  //          w   = w_f / alpha_here;                          // Return w with the conversion factor
  //          std::cout << "Our rounded w is: " << w << std::endl;
  //          grad_sample_points = loops_per_block * ceil(M_PI / (w * delta_t));
  //          grad_sample_points = ceil(grad_sample_points);
  //          points_per_loop    = (grad_sample_points / spectral_points) * 10; // check if it divi
  //      }

  //      // Recalculate omega (w) and sw based on above
  //      w  = M_PI * loops_per_block / (delta_t * grad_sample_points);
  //      sw = w / M_PI;

		rMrProt.getData().getsWipMemBlock().getalFree() [0] = static_cast<uint32_t>(points_per_loop); //Make available in the MDH
		rMrProt.getData().getsWipMemBlock().getalFree() [1] = static_cast<uint32_t>(ADCdwell); //Make available in the MDH
		
		 cout<<"The value of Nx is: "<<Nx<<endl;
		 cout<<"The value of Nsh is: "<<Nsh<<endl;
		//cout<<"The duration of water suppression pulse is:"<<wsatpulsedur<<endl;
		cout<<"The rf bandwidthh for excitation pulse (water) is: "<<rMrProt.getsSpecPara().getlRFExcitationBandwidth()<<endl;
		 cout<<"The value of fov is "<<fov<<endl;
         cout<<"The value of kmax is "<<kmax<<endl;
		 cout<<"The loops per block: "<<loops_per_block<<endl;
		 cout<<"The value of delta_t is "<<delta_t<<endl;
		cout<<"The number of points per loop is: "<<points_per_loop<<endl;
		 cout<<"The rounded value of gradient sample points is "<<grad_sample_points<<endl;
		cout<<"The original value of w is:"<<M_PI*sw<<endl;
		 cout<<"The rounded value of w is "<<w<<endl;
		 cout<<"The rounded value of spectral bandwidth is "<<sw<<endl;
		//cout<<"The original dt_adc is "<<dt_adc_max<<endl;
	

		// Prepping rosette gradients
		long gradsamplepoints= grad_sample_points;
		float *Gx= new float [gradsamplepoints];
		float *Gy= new float [gradsamplepoints];
		float *t=  new float [gradsamplepoints];
		cout<<"The value of w is: "<<w<<endl;

		for( int i=0; i<gradsamplepoints;i++)
		{ t[i]=static_cast<float>(i*delta_t);
		  Gx[i]= static_cast<float>(cos(2*w*t[i]));
		  Gy[i]= static_cast<float>(sin(2*w*t[i]));
	      
		} 
       
		float G_amp= static_cast<float>((w/Y)*kmax*1000); //mT/m
		cout<<"The value of G_amp is"<<G_amp;
		cout<<"The value of G_amp*Gx[0] is"<<G_amp*Gx[0];
		
	    // set ramp shape: array, ramp up points, ramp down points, last value NEZero: true/false
        m_RosGx.setRampShape(Gx, grad_sample_points, 0, true);
		m_RosGy.setRampShape(Gy, grad_sample_points,0, true);
		m_RosGx.setAmplitude(G_amp);
		m_RosGy.setAmplitude(G_amp); // should be multiplied by 1000 because sim shows mT/m
		m_RosGx.setDuration(10*grad_sample_points);
		m_RosGy.setDuration(10*grad_sample_points);
        

		 if( !m_RosGx.prep() )  {
          if( !(rSeqLim.isContextPrepForBinarySearch()) )
        {
			std::cout << "m_RosGx.prepExternal failed" << std::endl;
		}
              
        return m_RosGx.getNLSStatus();
        }
       

         if( !m_RosGy.prep() )  {
               if( !(rSeqLim.isContextPrepForBinarySearch()) )
        {
			std::cout << "m_RosGy.prepExternal failed" << std::endl;
		}
              return (m_RosGy.getNLSStatus());
        }


     //-----------------------------------------------------------------------
	 // Prepare Rosette Ramp up Gradients
	//----------------------------------------------------------------------
		
	
		int points_ramp=20;
		float Gx_rampup[20],Gy_rampup[20]; 
	
		for( int i=0; i<points_ramp;i++)
		{ 
		  Gx_rampup[i]= (i*Gx[0])/20 ;
		  Gy_rampup[i]= (i*Gy[0])/20 ;
	  
		} 
       
	    // set ramp shape: array, ramp up points, ramp down points, last value NEZero: true/false
        m_RosGx_rampup.setRampShape(Gx_rampup,points_ramp, 0, true);
		m_RosGy_rampup.setRampShape(Gy_rampup, points_ramp,0, true);
		m_RosGx_rampup.setAmplitude(G_amp);
		m_RosGy_rampup.setAmplitude(G_amp);
		m_RosGx_rampup.setDuration(10*points_ramp);
		m_RosGy_rampup.setDuration(10*points_ramp);
        

		if( !m_RosGx_rampup.prep() )  {
               if( !(rSeqLim.isContextPrepForBinarySearch()) )
        {
			std::cout << "m_RosGx_rampup.prepExternal failed" << std::endl;
		}
              return (m_RosGx_rampup.getNLSStatus());
        }
        
		if( !m_RosGy_rampup.prep() )  {
               if( !(rSeqLim.isContextPrepForBinarySearch()) )
        {
			std::cout << "m_RosGy_rampup.prepExternal failed" << std::endl;
		}
              return (m_RosGy_rampup.getNLSStatus());
        }
      
		
     //-----------------------------------------------------------------------
	 // Prepare Rosette Ramp down Gradients
	//----------------------------------------------------------------------
	    float *Gx_rampdown= new float [20];
		float *Gy_rampdown= new float [20];

		for (int j=0; j<points_ramp;j++)
		{ 
			Gx_rampdown[j]= Gx_rampup[points_ramp-j-1];
			Gy_rampdown[j]= Gy_rampup[points_ramp-j-1];
		}
        
        m_RosGx_rampdown.setRampShape(Gx_rampdown,points_ramp, 0, true);
		m_RosGy_rampdown.setRampShape(Gy_rampdown, points_ramp,0, true);
		m_RosGx_rampdown.setAmplitude(G_amp);
		m_RosGy_rampdown.setAmplitude(G_amp);
		m_RosGx_rampdown.setDuration(10*points_ramp);
		m_RosGy_rampdown.setDuration(10*points_ramp);
        	 
       	if( !m_RosGx_rampdown.prep() )  {
               if( !(rSeqLim.isContextPrepForBinarySearch()) )
        {
			std::cout << "m_RosGx_rampdown.prepExternal failed" << std::endl;
		}
              return (m_RosGx_rampdown.getNLSStatus());
        }
        
		if( !m_RosGy_rampdown.prep() )  {
               if( !(rSeqLim.isContextPrepForBinarySearch()) )
        {
			std::cout << "m_RosGy_rampdown.prepExternal failed" << std::endl;
		}
              return (m_RosGy_rampdown.getNLSStatus());
        }


    // gradient during excitation
    double dGradAmpl = 0.0;
    if ((dGradAmpl = m_rf_exc.getGSAmplitude()) > dGradMaxAmpl)
    {
        SEQ_TRACE_WARN.print("slice select. gradient of %f mT/m cannot be realized;", dGradAmpl);
        return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;
    }

    if (!SpecUtils::prepGradient(m_grad_exc, dGradAmpl, ((l_ExcDuration / 1000) + 1) * 1000, lRampTime_us)
        || !m_grad_exc.check())
        return m_grad_exc.getNLSStatus();

    // the refocussing gradient
    // is balanced with the spoiler after the 1st refocussing pulse

    // gradient during refocusing pulse
    // the refocusing gradient is balanced with the spoiler after the 1st refocussing pulse

    // if(  (dGradAmpl = m_rf_ref.getGSAmplitude()) > dGradMaxAmpl )
    //     {
    //         SEQ_TRACE_WARN.print("slice select. gradient of %f mT/m cannot be realized;", dGradAmpl);
    //         return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;
    //     }
    

    // if (!SpecUtils::prepGradient(m_grad_ref, dGradAmpl, l_RefDuration , lRampTime_us)
    //     || !m_grad_ref.check())
    //     return m_grad_ref.getNLSStatus();

    // compute phase encoding gradients
    // 1st phase encoding direction is READOUT
    m_d_1st_csi_grad_offset = 0; // no offset due to spoiling etc.
    m_d_1st_csi_grad_step
        = SpecUtils::calcCSIGradStep(dLarmorConst, rMrProt.sliceSeries().front().readoutFOV(), encoddur);

    // do some checking
    dmin
        = m_d_1st_csi_grad_offset - m_d_1st_csi_grad_step * std::lround(0.5 * rMrProt.kSpace().baseResolution());
    dmax
        = m_d_1st_csi_grad_offset + m_d_1st_csi_grad_step * std::lround(0.5 * rMrProt.kSpace().baseResolution());

    // gradient overflow
    if (fabs(dmin) > dGradMaxAmpl || fabs(dmax) > dGradMaxAmpl)
    {
        SEQ_TRACE_WARN.print(
            "phase encod. gradient READOUT dir. of %f mT/m cannot be realized;",
            (fabs(dmin) > fabs(dmax) ? fabs(dmin) : fabs(dmax)));
        return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;
    }

    // 2nd phase encoding direction is PHASE
    m_d_2nd_csi_grad_offset = 0; // no offset due to spoiling etc.
    m_d_2nd_csi_grad_step
        = SpecUtils::calcCSIGradStep(dLarmorConst, rMrProt.sliceSeries().front().phaseFOV(), encoddur);

    // do some checking
    dmin = m_d_2nd_csi_grad_offset - m_d_2nd_csi_grad_step * std::lround(0.5 * rMrProt.kSpace().phaseEncodingLines());
    dmax = m_d_2nd_csi_grad_offset + m_d_2nd_csi_grad_step * std::lround(0.5 * rMrProt.kSpace().phaseEncodingLines());

    // gradient overflow
    if (fabs(dmin) > dGradMaxAmpl || fabs(dmax) > dGradMaxAmpl)
    {
        SEQ_TRACE_WARN.print(
            "phase encod. gradient PHASE dir. of %f mT/m cannot be realized;",
            (fabs(dmin) > fabs(dmax) ? fabs(dmin) : fabs(dmax)));
        return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;
    }

    // 3rd phase encoding direction is SLICE
    m_d_3rd_csi_grad_offset = -(refoc * m_rf_exc.getDuration() + .5 * m_grad_exc.getRampDownTime())
                              * m_rf_exc.getGSAmplitude()
                              / (double)encoddur; // gradient strength offset due to slice rephasing

    m_asymAmp = - (refoc * m_rf_exc.getDuration() + 0.5 * m_grad_exc.getRampDownTime() ) * m_rf_exc.getGSAmplitude() / (double)encoddur;
                              
    m_d_3rd_csi_grad_step
        = SpecUtils::calcCSIGradStep(dLarmorConst, rMrProt.sliceSeries().aFront().thickness(), encoddur);

    // do some checking
    // the phase encoding gradients in the 3rd direction are only relevant if we perform a 3D measurement
    if (rMrProt.kSpace().getlPartitions() > 1)
    {
        dmin = m_d_3rd_csi_grad_offset - m_d_3rd_csi_grad_step * std::lround(0.5 * rMrProt.kSpace().getlPartitions());
        dmax = m_d_3rd_csi_grad_offset + m_d_3rd_csi_grad_step * std::lround(0.5 * rMrProt.kSpace().getlPartitions());
    }
    else
    {
        dmin = m_d_3rd_csi_grad_offset;
        dmax = m_d_3rd_csi_grad_offset;
    }

    // gradient overflow
    if (fabs(dmin) > dGradMaxAmpl || fabs(dmax) > dGradMaxAmpl)
    {
        SEQ_TRACE_ALWAYS.print(
            "phase encod. gradient SLICE dir. of %f mT/m cannot be realized;",
            (fabs(dmin) > fabs(dmax) ? fabs(dmin) : fabs(dmax)));
        return MRI_SEQ_SEQU_SEQ_NOT_PREPARED;
    }
    // gradient ramping
    // before: slice encoding, after: nullptr
    // BUT: since the slice select. gradient and the refoc. gradient each have their own ramp durations,
    // the grad_max checking above is sufficient

    // prepare encoding gradient timing
    // the amplitude is set within fSeqRun()
    m_encod_sl.set(lRampTime_us, encoddur, lRampTime_us);
    // m_encod_ph.set(lRampTime_us, encoddur, lRampTime_us);
    // m_encod_ro.set(lRampTime_us, encoddur, lRampTime_us);

    // calculate the number of requests of each kernel call
    long lKernelRequestsPerMeasurement = rMrProt.getsSpecPara().getlPreparingScans();

    // check dimensions
    if (rMrProt.getsSpecPara().getlFinalMatrixSizeRead() < rMrProt.kSpace().getlBaseResolution()
        || rMrProt.getsSpecPara().getlFinalMatrixSizePhase() < rMrProt.kSpace().getlPhaseEncodingLines()
        || ((rMrProt.kSpace().dimension() == SEQ::DIM_3) ? rMrProt.getsSpecPara().getlFinalMatrixSizeSlice() : 1)
               < rMrProt.kSpace().getlPartitions())
        return MRI_SEQ_SEQU_ERROR;

    // compute the CSI gradient table
    SpecUtils::calcCSIGradientTable(
        rMrProt.kSpace().getlBaseResolution(),
        rMrProt.kSpace().getlPhaseEncodingLines(),
        rMrProt.kSpace().getlPartitions(),
        rMrProt.averages(),
        rMrProt.getsSpecPara().getlPhaseEncodingType() == SEQ::PHASE_ENCODING_FULL,
        rMrProt.getsSpecPara().getlPhaseEncodingType() == SEQ::PHASE_ENCODING_WEIGHTED,
        m_sh_1st_csi_addr,
        m_sh_2nd_csi_addr,
        m_sh_3rd_csi_addr,
        m_sh_csi_weight,
        lKernelRequestsPerMeasurement);

    lKernelRequestsPerMeasurement = rMrProt.getsSpecPara().getlPreparingScans(); //repeat this line to have correct meas time (temporary solution, EK)

    if (!(m_encod_sl.prepAmplitude( m_asymAmp ))
        || !(m_encod_sl.check()))
        return m_encod_sl.getNLSStatus();

    // if (!(m_encod_ro.prepAmplitude(m_d_1st_csi_grad_offset + m_sh_1st_csi_addr[0] * m_d_1st_csi_grad_step))
    //     || !(m_encod_ro.check()))
    //     return m_encod_ro.getNLSStatus();

    // if (!(m_encod_ph.prepAmplitude(m_d_2nd_csi_grad_offset + m_sh_2nd_csi_addr[0] * m_d_2nd_csi_grad_step))
    //     || !(m_encod_ph.check()))
    //     return m_encod_ph.getNLSStatus();

//////////////////////////////////////////////////////////////////////////
	// Prephaser gradients
	//////////////////////////////////////////////////////////////////////////
	    
	// Prephaser gradients	
	double PrephaseMoment_x = 0.5*m_RosGx_rampup.getTotalTime()*m_RosGx_rampup.getAmplitude();
	cout<<"The amplitude of ramp up gradient is"<<m_RosGx_rampup.getAmplitude()<<endl;
	cout<<"The max magnitude of ramp up gradient "<<m_RosGx_rampup.getMaxMagnitude()<<endl;

    m_RosGx_Prephase.setRampUpTime( m_encod_sl.getRampUpTime());
    m_RosGx_Prephase.setRampDownTime( m_encod_sl.getRampDownTime());
    m_RosGx_Prephase.setDuration( m_encod_sl.getDuration()); 
    m_RosGx_Prephase.setAmplitude(-PrephaseMoment_x/m_RosGx_Prephase.getDuration());

    if( !m_RosGx_Prephase.prep() )  {
                if( !(rSeqLim.isContextPrepForBinarySearch()) )
        {
            std::cout << "m_RosGx_Prephase.prepExternal failed" << std::endl;
        }
            return (m_RosGy_Prephase.getNLSStatus());
        }    

        double PrephaseMoment_y = 0.5*m_RosGy_rampup.getTotalTime()*0*m_RosGy_rampup.getAmplitude();
	
        m_RosGy_Prephase.setRampUpTime( m_encod_sl.getRampUpTime());
        m_RosGy_Prephase.setRampDownTime( m_encod_sl.getRampDownTime());
        m_RosGy_Prephase.setDuration( m_encod_sl.getDuration()); 
        m_RosGy_Prephase.setAmplitude(-PrephaseMoment_y/m_RosGy_Prephase.getDuration());
   
        if( !m_RosGy_Prephase.prep() )  {
              if( !(rSeqLim.isContextPrepForBinarySearch()) )
             {
               std::cout << "m_RosGy_Prephase.prepExternal failed" << std::endl;
           }
              return (m_RosGy_Prephase.getNLSStatus());
           }  

    // // spoiler 1

    // //m_sp1_ph.setMaxMagnitude(max_grad_ampl);
    // if( !( m_sp1_ph.prepAmplitude( lRampTime_us, sp1dur, lRampTime_us, sp1ampl ) ) || !( m_sp1_ph.check() ) )
    // {
    //     std::cout << "m_sp1_ph.prepAmplitude failed" << std::endl;
    //     return m_sp1_ph.getNLSStatus();
    // }

    // //m_sp1_ro.setMaxMagnitude(max_grad_ampl);
    // if( !( m_sp1_ro.prepAmplitude( lRampTime_us, sp1dur, lRampTime_us, sp1ampl ) ) || !( m_sp1_ro.check() ) )
    // {
    //     std::cout << "m_sp1_ro.prepAmplitude failed" << std::endl;
    //     return m_sp1_ro.getNLSStatus();
    // }

    // //m_sp1_sl.setMaxMagnitude(max_grad_ampl);
    // if( !( m_sp1_sl.prepAmplitude( lRampTime_us, sp1dur, lRampTime_us, ((sp1ampl*sp1dur)+(m_encod_sl.getAmplitude()*m_encod_sl.getDuration()))/sp1dur ) ) || !( m_sp1_sl.check() ) )
    // {
    //     std::cout << "m_sp1_sl.prepAmplitude failed" << std::endl;
    //     return m_sp1_sl.getNLSStatus();
    // }

    // // spoiler 2

    // //m_sp2_ph.setMaxMagnitude(max_grad_ampl);
    // if( !( m_sp2_ph.prepAmplitude( lRampTime_us, sp2dur, lRampTime_us, ((sp2ampl*sp2dur)+(m_RosGy_Prephase.getDuration()*m_RosGy_Prephase.getAmplitude()))  / sp2dur ) ) || !( m_sp2_ph.check() ) )
    // {
    //        std::cout << "m_sp2_ph.prepAmplitude failed" << std::endl;
    //       return m_sp2_ph.getNLSStatus();
    // }

    // //m_sp2_ro.setMaxMagnitude(max_grad_ampl);
    // if( !( m_sp2_ro.prepAmplitude( lRampTime_us, sp2dur, lRampTime_us, ((sp2ampl*sp2dur)+(m_RosGx_Prephase.getDuration()*m_RosGx_Prephase.getAmplitude()))  / sp2dur ) ) || !( m_sp2_ro.check() ) )
    // {
    //        std::cout << "m_sp2_ro.prepAmplitude failed" << std::endl;
    //       return m_sp2_ro.getNLSStatus();
    // }

    // //m_sp2_sl.setMaxMagnitude(max_grad_ampl);
    // if( !( m_sp2_sl.prepAmplitude( lRampTime_us, sp2dur, lRampTime_us, sp2ampl ) ) || !( m_sp2_sl.check() ) )
    // {
    //        std::cout << "m_sp2_sl.prepAmplitude failed" << std::endl;
    //       return m_sp2_sl.getNLSStatus();
    // }


	  // --------------------------------------------------------------------
      // Calculate TEFill-times and check, whether timing can be realized
      // --------------------------------------------------------------------
     //lb template for delays: m_lTEFill_us =rMrProt.te()[0] - static_cast<long>(m_rf_exc.getDuration() *(1- m_rf_exc.getAsymmetry())) - m_encod_sl.getDuration()- m_encod_sl.getRampDownTime()- static_cast<long>(m_RosGx_rampup.getTotalTime())- static_cast<long>(m_grad_exc.getRampDownTime());
     m_lTEFill_us =rMrProt.te()[0] - static_cast<long>(m_rf_exc.getDuration() *(1- m_rf_exc.getAsymmetry())) - m_encod_sl.getDuration()- m_encod_sl.getRampDownTime()- static_cast<long>(m_RosGx_rampup.getTotalTime())- static_cast<long>(m_grad_exc.getRampDownTime());

    // std::cout<<"m_rf_exc.getDuration() *(1- m_rf_exc.getAsymmetry())"<<static_cast<long>(m_rf_exc.getDuration() *(1- m_rf_exc.getAsymmetry()))<<std::endl;
    // std::cout<<"m_encod_sl.getDuration()"<<m_encod_sl.getDuration()<<std::endl;
    // std::cout<<"m_encod_sl.getRampDownTime()"<<m_encod_sl.getRampDownTime()<<std::endl;
    // std::cout<<" static_cast<long>(m_RosGx_rampup.getTotalTime())"<< static_cast<long>(m_RosGx_rampup.getTotalTime())<<std::endl;
    // std::cout<<"static_cast<long>(m_grad_exc.getRampDownTime())"<<static_cast<long>(m_grad_exc.getRampDownTime())<<std::endl;
    // std::cout<<"m_lTEFill_us"<<m_lTEFill_us<<std::endl;
    
	 //delay_1 = rMrProt.te()[0]/2 - m_rf_exc.getDuration()/2 -m_grad_exc.getRampDownTime() - m_encod_sl.getTotalTime() - m_grad_ref.getRampUpTime() - m_rf_ref.getDuration()/2;
	//  delay_1 = rMrProt.te()[0]/2 - m_rf_exc.getDuration()/2 -m_grad_exc.getRampDownTime() - m_sp1_sl.getTotalTime() - (m_grad_ref.getDuration() - m_rf_ref.getDuration() ) - m_rf_ref.getDuration()/2;
    //  delay_1 = rMrProt.te()[0]/2 - m_rf_exc.getDuration()/2 -m_grad_exc.getRampDownTime() /*- m_sp1_sl.getTotalTime()*/;
	//  delay_1 = rMrProt.te()[0]/2 - (m_rf_exc.getDuration()/2) - m_sp1_sl.getDuration() - m_grad_ref.getRampUpTime() - (m_rf_ref.getDuration()/2);
	//  std::cout<<"The value of TE is inputted as: "<<rMrProt.te()[0]<<std::endl;
	//  std::cout<<"The value of m_rf_exc.getDuration()/2 is inputted as: "<<m_rf_exc.getDuration()/2<<std::endl;
	//  std::cout<<"The value of m_grad_exc.getRampDownTime() is inputted as: "<<m_grad_exc.getRampDownTime()<<std::endl;
	//  std::cout<<"The value of m_encod_sl.getTotalTime() is inputted as: "<<m_encod_sl.getTotalTime()<<std::endl;
	//  std::cout<<"The value of m_grad_ref.getRampUpTime() is inputted as: "<<m_grad_ref.getRampUpTime()<<std::endl;
	//  std::cout<<"The value of m_rf_ref.getDuration() is inputted as: "<<m_rf_ref.getDuration()/2<<std::endl;
	//  std::cout<<"The value of delay_1 is: "<<delay_1<<std::endl;

	//  if ( delay_1<0 )
    if (m_lTEFill_us<0)
	 { 
		  if( !(rSeqLim.isContextPrepForBinarySearch()) )
          {
		    std::cout << "TE value not possible. Aborting!!" << std::endl;
		}
           return MRI_SEQ_SEQU_ERROR;
	 }

	 //delay_2 = rMrProt.te()[0]/2 - m_rf_ref.getDuration()/2 - m_grad_ref.getRampDownTime() - m_sp2_sl.getTotalTime() - m_RosGx_rampup.getDuration();
	//  delay_2 = rMrProt.te()[0]/2 - (m_rf_ref.getDuration()/2) - m_sp2_sl.getDuration() - m_RosGx_rampup.getDuration();
	//  std::cout<<"The value of delay_2 is: "<<delay_2<<std::endl;
	//  if ( delay_2<0 )
	//  { 
	// 	  if( !(rSeqLim.isContextPrepForBinarySearch()) )
    //       {
	// 	    std::cout << "TE value not possible. Aborting!!" << std::endl;
	// 	}
    //        return MRI_SEQ_SEQU_ERROR;
	//  }


    delay_3 = ceil((2*M_PI) / (n_ti*w*2)*1000000);
    std::cout << "The value of delay_3 is: " << delay_3 << std::endl;

    if (delay_3 < 0)
    {
        if (!(rSeqLim.isContextPrepForBinarySearch()))
        {
            std::cout << "Delay_3 value not possible. Aborting!!" << std::endl;
        }
        return MRI_SEQ_SEQU_ERROR;
    }

	 //cout<<"Hence I calculate TE = "<< m_rf_exc.getDuration()/2 + m_grad_exc.getRampDownTime() + m_encod_sl.getDuration() + m_enco<<endl;
	 /////////////////////////////////////////////////////////////////////////
	 //Calculation of new KernelRequestsPerMeasurement
	 ////////////////////////////////////////////////////////////////////////
	 lKernelRequestsPerMeasurement += static_cast<long>(Nsh*(rMrProt.averages())*n_ti);
   
    


    // Calculate the total measurement time, including measurement repeats
    const double dMeasureTimeUsec = static_cast<double>(lKernelRequestsPerMeasurement) * rMrProt.tr()[0];

    double dTotalMeasureTimeMsec = 0.0;
    lStatus = fSBBMeasRepetDelaysPrep(rMrProt, rSeqLim, rSeqExpo, dMeasureTimeUsec / 1000.0, &dTotalMeasureTimeMsec);
    CheckStatusPR(lStatus, "fSBBMeasRepetDelaysPrep");

    // Output DMeas variables
	cout<<" The value of dMeasureTimeUsec is "<<dMeasureTimeUsec<<endl;
	cout<<" The value of dTotalMeasureTimeMsec is "<<dTotalMeasureTimeMsec<<endl;
	cout<<" The value of lKernelRequestsPerMeasurement is: "<<lKernelRequestsPerMeasurement<<endl;
	cout<<" The value of the current TR is :"<<rMrProt.tr()[0]<<endl;


    // include the repetitions to the number of kernel requests AFTER calculating dTotalMeasureTimeMsec,
    // since otherwise the repetitions would have been taken into account twice
    const long lKernelRequestsTotal = lKernelRequestsPerMeasurement * (rMrProt.repetitions() + 1);
    
    // final spoiling pulses
    const long tau = 20000;

    if (!(m_finsp_ph.prepAmplitude(lRampTime_us, tau, lRampTime_us, 5.0)) || !(m_finsp_ph.check()))
        return m_finsp_ph.getNLSStatus();

    if (!(m_finsp_ro.prepAmplitude(lRampTime_us, tau, lRampTime_us, 5.0)) || !(m_finsp_ro.check()))
        return m_finsp_ro.getNLSStatus();

    if (!(m_finsp_sl.prepAmplitude(lRampTime_us, tau, lRampTime_us, 5.0)) || !(m_finsp_sl.check()))
        return m_finsp_sl.getNLSStatus();

    const long lFinalSpoilDuration = tau + lRampTime_us;

    // Prepare Sequence Building block (SBBWatSat)
    if (!m_sWaterSat.prep(rMrProt, rSeqLim, rSeqExpo))
        return m_sWaterSat.getNLSStatus(); // only if main nucleus = 1H

    const long lSBBDuration = m_sWaterSat.getDurationPerRequest();

    rfInfoMainNucPerScan += m_sWaterSat.getRFInfoPerRequest();

    // Calculate various start times
    // m_adc1.setStartTime(m_grad_exc.getTotalTime() + m_encod_sl.getTotalTime()); // encoding duration defines TE
    // if (rMrProt.getsSpecPara().getlNOEType() == SEQ::NOE_RECTANGULAR)
    // {
    //     m_adc1.setStartTime(m_adc1->getStartTime() + noe_dur);
    // }

    // Calculate minimum TR
    // const long lMinTR = fSDSRoundUpGRT(
    //     SysProperties::getCoilCtrlLead() + lSBBDuration + m_adc1.getStartTime() + m_adc1.getRoundedDuration()
    //     + m_finsp_ph.getTotalTime());

    // m_lTRFill = rMrProt.tr()[0] - lMinTR;
    // if (m_lTRFill < 0)
    //     return MRI_SEQ_SEQU_NEGATIV_TRFILL;

    // set the receiver gain
    lStatus = m_rxGainSelector.setRXGainCodeToSystem(rMrProt, rSeqLim);
    CheckStatusPR(lStatus, "setRXGainCodeToSystem");

    // for estimation of residual measurement time
    rMrProt.physiology().getPhysioMode(m_firstSignal, m_firstMethod, m_secondSignal, m_secondMethod);

    long lPhysioHalts = 0;
    if (rSeqLim.isContextPrepForScanTimeCalculation() || rSeqLim.isContextNormal())
    {
        if (m_firstMethod == SEQ::METHOD_TRIGGERING)
        {
            lPhysioHalts = rMrProt.averages() * static_cast<long>(m_sh_1st_csi_addr.size());
            lPhysioHalts *= (rMrProt.repetitions() + 1);
        }
        else
        {
            lPhysioHalts = 0;
        }
        SeqUT.setExpectedPhysioHalts(lPhysioHalts);
    }

    if (m_firstMethod == SEQ::METHOD_TRIGGERING)
    {
        m_lKernelCallsPerRelevantSignal
            = std::max(1L, static_cast<long>(lKernelRequestsPerMeasurement * 1000000.0 / dMeasureTimeUsec)); // per sec
        rSeqExpo.setRelevantReadoutsForMeasTime(
            (lKernelRequestsPerMeasurement - rMrProt.getsSpecPara().getlPreparingScans())
            / m_lKernelCallsPerRelevantSignal);
    }

    // checking of sequence & output of MRI_SEQ_SEQU_ERROR to actuate solve handler
    // VecSizeTRConflict & BandWidthTRConflict
    m_lTRNeededSpectro = n_blocks*m_adc1[0].getRoundedDuration() + lSBBDuration + lFinalSpoilDuration + m_grad_exc.getTotalTime() + m_encod_sl.getTotalTime()
                         /*+ m_sp1_sl.getTotalTime() + m_grad_ref.getTotalTime()*/ + m_RosGx_rampup.getTotalTime() + noe_dur + 1000;

    if (m_lTRNeededSpectro > rMrProt.tr()[0])
    {
        if (!(rSeqLim.isContextPrepForBinarySearch()))
        {
            SEQ_TRACE_WARN.print("TR = %d us cannot be realized;", rMrProt.tr()[0]);
        }
        return MRI_SEQ_SEQU_ERROR;
    }

    // in case of PhysioTrigger add PulseDuration to MeasTime
    if (lPhysioHalts > 0)
    {
        long l1           = (long)(rMrProt.physiology().triggerDelay(m_firstSignal));
        long l2           = getMinSBBECGFillTimeRunDelay(*m_firstHalt); // retrieve min duration of SBBECGFillTimeRun
        long triggerDelay = l1 > l2 ? l1 : l2;

        dTotalMeasureTimeMsec += static_cast<double>(lPhysioHalts * triggerDelay) / 1.e3;
    }

    // export parameters to ICE program
    rSeqExpo.resetAllRFBlockInfos();

    if (txNucleus.isValid())
    {
        rSeqExpo.addRFBlockInfo(
            dMeasureTimeUsec * 1.0e-6, // [us] -> [s]
            SeqExpoRFBlockInfo::VALUETYPE_ACTUAL,
            mainNucleus.get(),
            lKernelRequestsTotal * rfInfoMainNucPerScan,
            SeqExpoRFBlockInfo::VALUETYPE_ACTUAL,
            txNucleus.get(),
            lKernelRequestsTotal * rfInfoTxNucPerScan,
            SeqExpoRFBlockInfo::VALUETYPE_ACTUAL);
    }
    else
    {
        rSeqExpo.addRFBlockInfo(
            dMeasureTimeUsec * 1.0e-6, // [us] -> [s]
            SeqExpoRFBlockInfo::VALUETYPE_ACTUAL,
            mainNucleus.get(),
            lKernelRequestsTotal * rfInfoMainNucPerScan,
            SeqExpoRFBlockInfo::VALUETYPE_ACTUAL);
    }

    rSeqExpo.setMeasureTimeMin(dMeasureTimeUsec / 60000000.0);
    rSeqExpo.setTotalMeasureTimeMin(dTotalMeasureTimeMsec / 60000.0);

    // rSeqExpo.setMeasuredPELines( 1 );
    rSeqExpo.setSequenceString("csi_fid");
    rSeqExpo.setSeqShortString("csi_fid");
    rSeqExpo.setMeasuredPELines(rMrProt.kSpace().phaseEncodingLines()); // does not account for elliptical or weighted
    rSeqExpo.setMeasured3dPartitions(rMrProt.kSpace().partitions());    // does not account for elliptical or weighted

    return lStatus;
}

NLSStatus Csi_fid::check(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, SEQCheckMode*)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    // execute kernel for checking (GSWD look ahead functionality)
    if (!(m_encod_sl.prepAmplitude( m_asymAmp ))
        || !(m_encod_sl.check()))
        return m_encod_sl.getNLSStatus();

    // if (!(m_encod_ro.prepAmplitude(m_d_1st_csi_grad_offset + m_sh_1st_csi_addr[0] * m_d_1st_csi_grad_step))
    //     || !(m_encod_ro.check()))
    //     return m_encod_ro.getNLSStatus();

    // if (!(m_encod_ph.prepAmplitude(m_d_2nd_csi_grad_offset + m_sh_2nd_csi_addr[0] * m_d_2nd_csi_grad_step))
    //     || !(m_encod_ph.check()))
    //     return m_encod_ph.getNLSStatus();

    lStatus = runKernel(rMrProt, rSeqLim, rSeqExpo, KERNEL_CHECK);
    CheckStatusPR(lStatus, "runKernel");

    lStatus = runKernel(rMrProt, rSeqLim, rSeqExpo, KERNEL_CHECK);
    CheckStatusPR(lStatus, "runKernel");

    return MRI_SEQ_SEQU_NORMAL;
}

NLSStatus Csi_fid::run(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    int Nx= rMrProt.kSpace().baseResolution();
	int Nsh=  static_cast<int>(ceil((M_PI*Nx)/2));

    //  initialization of the unit test function
    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunStart, 0, 0, 0, 0, 0);

    //  set looping parameters
    // uint16_t d1 = static_cast<uint16_t>(rMrProt.kSpace().baseResolution()) / 2;
    // uint16_t d2 = static_cast<uint16_t>(rMrProt.kSpace().phaseEncodingLines()) / 2;
    // uint16_t d3 = static_cast<uint16_t>(rMrProt.kSpace().partitions()) / 2;

    for(int m=0;m<9;m++)
	    {
        m_adc1[m].getMDH().setCseg(0); // 1st
        m_adc1[m].getMDH().setClin(0); // 2nd
        m_adc1[m].getMDH().setCpar(0); // 3rd
        m_adc1[m].getMDH().setCphs(0);
        m_adc1[m].getMDH().setCeco(0); // echo number
        m_adc1[m].getMDH().setCslc(0); // slice number
        m_adc1[m].getMDH().setCrep(0); // repetitions

        // other Mdh info
        m_adc1[m].getMDH().setIceProgramPara(4, static_cast<uint16_t>(m_lSamplesBeforeEcho));
        }
    // execute repetition loop
    long       lCurrKernelCalls = 0;
    const long n_meas           = rMrProt.measurements();
    for (long k = 0; k < n_meas; k++)
    {
        for(int m=0;m<9;m++)
		{
        m_adc1[m].getMDH().setCrep(static_cast<uint16_t>(k));
        }

        // execute prepare loop
        fRTSetReadoutEnable(0); // disable ADC events
        for (long i = 0; i < rMrProt.getsSpecPara().getlPreparingScans(); i++)
        {
            lStatus = runKernel(rMrProt, rSeqLim, rSeqExpo, KERNEL_IMAGE);
            CheckStatusPR(lStatus, "runKernel");
        }
        fRTSetReadoutEnable(1); // enable ADC events

        // execute acquisition loop
        // long lCurrKernelCalls = 0;

        SEQ_TRACE_WARN.print("HEHE00");
        if (rMrProt.kSpace().averagingMode() == SEQ::INNER_LOOP) // acquisition inside lines loop / Short Term / normal
        {
            // nave = rMrProt.averages();
            for (size_t i = 0; i < static_cast<size_t>(Nsh); i++)
            {
                if (!(m_encod_sl.prepAmplitude( m_asymAmp ))
                    || !(m_encod_sl.check()))
                    return m_encod_sl.getNLSStatus();

                // if (!(m_encod_ro.prepAmplitude(m_d_1st_csi_grad_offset + m_sh_1st_csi_addr[i] * m_d_1st_csi_grad_step))
                //     || !(m_encod_ro.check()))
                //     return m_encod_ro.getNLSStatus();

                // if (!(m_encod_ph.prepAmplitude(m_d_2nd_csi_grad_offset + m_sh_2nd_csi_addr[i] * m_d_2nd_csi_grad_step))
                //     || !(m_encod_ph.check()))
                //     return m_encod_ph.getNLSStatus();

                SEQ_TRACE_WARN.print("HEHE001");
                for (long j = 0; j < rMrProt.averages(); j++) // averages
                {
                    for (long l = 0; l < n_ti; l++)
                    {
                        SEQ_TRACE_WARN.print("HEHE0002");
                        m_ISISalternator = l % n_ti;
                        SEQ_TRACE_WARN.print("HEHE003");
                        for (int m = 0; m < 9; m++)
                        {
                            m_adc1[m].getMDH().setClin(static_cast<uint16_t>(i));
                        }
                        sNewRotationMatrix = m_fov->getROT_MATRIX();
                        dRotAngle          = ((2.0 * M_PI) / Nsh) * i;
                        ZRotMatrix[0][0]   = cos(dRotAngle);
                        ZRotMatrix[0][1]   = -sin(dRotAngle);
                        ZRotMatrix[0][2]   = 0;
                        ZRotMatrix[1][0]   = sin(dRotAngle);
                        ZRotMatrix[1][1]   = cos(dRotAngle);
                        ZRotMatrix[1][2]   = 0;
                        ZRotMatrix[2][0]   = 0;
                        ZRotMatrix[2][1]   = 0;
                        ZRotMatrix[2][2]   = 1;
                        sExistingMatrix    = m_fov->getROT_MATRIX();
                        // Matrix Multiplication
                        for (lJ = 0; lJ < 3; lJ++)
                        {
                            for (lK = 0; lK < 3; lK++)
                            {
                                sNewRotationMatrix.dMat[lJ][lK] = 0.0;
                                for (lL = 0; lL < 3; lL++)
                                {
                                    sNewRotationMatrix.dMat[lJ][lK]
                                        += sExistingMatrix.dMat[lJ][lL] * ZRotMatrix[lL][lK];
                                }
                            }
                        }
                        for (int m = 0; m < 9; m++)
                        {
                            m_adc1[m].getMDH().setCacq(static_cast<uint16_t>(j)); // averages
                            // m_adc1[m].getMDH().setCseg(static_cast<uint16_t>(
                            //     m_sh_1st_csi_addr[i] + d1)); // Cseg must be between 0 and
                            //     kSpace().phaseEncodingLines() - 1
                            // m_adc1[m].getMDH().setClin(static_cast<uint16_t>(
                            //     m_sh_2nd_csi_addr[i] + d2)); // Clin must be between 0 and
                            //     kSpace().phaseEncodingLines() - 1
                            // m_adc1[m].getMDH().setCpar(static_cast<uint16_t>(
                            //     m_sh_3rd_csi_addr[i] + d3)); // Cpar must be between 0 and kSpace().partitions() - 1
                            SEQ_TRACE_WARN.print("HEHE004");
                            m_adc1[m].getMDH().setCeco(l); // preparing temporal interleaves CF
                            SEQ_TRACE_WARN.print("HEHE005");
                            // flags for extracting time-stamps
                            m_adc1[m].getMDH().setFirstScanInSlice(!i && !j && !l);
                            m_adc1[m].getMDH().setLastScanInSlice(
                                i == (static_cast<size_t>(Nsh) - 1) && j == (rMrProt.averages() - 1) && l == (n_ti - 1));
                            m_adc1[m].getMDH().setLastScanInConcat(
                                i == (static_cast<size_t>(Nsh) - 1)
                                && j == (rMrProt.averages() - 1)); // flag for concatenations
                            m_adc1[m].getMDH().setLastScanInMeas(
                                i == (static_cast<size_t>(Nsh) - 1) && j == (rMrProt.averages() - 1) && l == (n_ti - 1)); // end of scan
                            SEQ_TRACE_WARN.print("HEHE006");
                        }

                        // realize triggering
                        rMrProt.physiology().getPhysioMode(
                            m_firstSignal, m_firstMethod, m_secondSignal, m_secondMethod);
                        if (m_firstMethod == SEQ::METHOD_TRIGGERING /* && m_eTriggerMode == SEQ::TRIGGER_STANDARD */)
                        {
                            long l1 = rMrProt.physiology().triggerDelay(m_firstSignal);
                            long l2 = getMinSBBECGFillTimeRunDelay(
                                *m_firstHalt); // retrieve min duration of SBBECGFillTimeRun
                            const long triggerDelay = l1 > l2 ? l1 : l2;
                            const long delayPriorToTrigger
                                = 0; // no gradient ramps must fall within the trigger event block

                            CheckStatusPR(
                                lStatus = fSBBECGFillTimeRun(&m_firstHalt, triggerDelay, delayPriorToTrigger),
                                "fSBBECGFillTimeRun");

                            // for estimation  of residual measurement time
                            lCurrKernelCalls++;
                            if (!(lCurrKernelCalls % m_lKernelCallsPerRelevantSignal))
                            {
                                for (int m = 0; m < 9; m++)
                                {
                                    m_adc1[m].setRelevantForMeasTime();
                                }
                            }
                        }

                        // realize off-centre FoV positions by incrementing the phase of the excitation pulses from step
                        // to step
                        double excit_phase = SpecUtils::calcCSIPhase(
                            m_sh_1st_csi_addr[i],
                            m_d_read_pos,
                            rMrProt.sliceSeries().front().readoutFOV(),
                            m_sh_2nd_csi_addr[i],
                            m_d_phase_pos,
                            rMrProt.sliceSeries().front().phaseFOV(),
                            m_sh_3rd_csi_addr[i],
                            m_d_slice_pos,
                            rMrProt.sliceSeries().front().thickness());

                        // double aqc_phase = 0.0;

                        if ((m_sh_1st_csi_addr[i] ^ m_sh_2nd_csi_addr[i] ^ m_sh_3rd_csi_addr[i] ^ j) & 1)
                        {
                            excit_phase += 180.0;
                            // aqc_phase = 180.0;
                        }
                        // else
                        // aqc_phase = 0.0;

                        // m_ph_s_exc.increasePhase(excit_phase);
                        // m_ph_n_exc.decreasePhase(excit_phase);
                        // m_ph_s_adc.increasePhase(aqc_phase);
                        // m_ph_n_adc.decreasePhase(aqc_phase);
                        SEQ_TRACE_WARN.print("HEHE008");
                        lStatus = runKernel(rMrProt, rSeqLim, rSeqExpo, KERNEL_IMAGE);
                        SEQ_TRACE_WARN.print("HEHE009");

                        
                        CheckStatusPR(lStatus, "runKernel");
                        SEQ_TRACE_WARN.print("HEHE010");

                        // undo phase cycling
                        // m_ph_s_exc.decreasePhase(excit_phase);
                        // m_ph_n_exc.increasePhase(excit_phase);
                        // m_ph_s_adc.decreasePhase(aqc_phase);
                        // m_ph_n_adc.increasePhase(aqc_phase);
                    }
                } // end averaging loop
            }     // end encoding loop
        } // if (rMrProt.kSpace().averagingMode() == SEQ::INNER_LOOP)    // acquisition inside lines loop / Short Term /
          // normal

        SEQ_TRACE_WARN.print("HEHE1");
        if (rMrProt.kSpace().averagingMode() == SEQ::OUTER_LOOP) // acquisition outside lines loop
        {
            for (long j = 0; j < rMrProt.averages(); j++) // averages
            {
                if (!(m_encod_sl.prepAmplitude(m_asymAmp)) || !(m_encod_sl.check()))
                    return m_encod_sl.getNLSStatus();

                for (size_t i = 0; i < static_cast<size_t>(Nsh); i++) // PE steps
                {
                    SEQ_TRACE_WARN.print("HEHE2");
                    for (long l = 0; l < n_ti; l++)
                    {
                        SEQ_TRACE_WARN.print("HEHE3");
                        m_ISISalternator = l % n_ti;
                        #ifdef _DVP_DEBUG
                        std::cout << std::endl;
                        // cout << "Kind of averaging" << endl;
                        std::cout << "Average: " << j + 1 << ";   j: " << j << std::endl;
                        // cout << "max. Average: nave:" << nave << endl;
                        // cout << "PE step:" << i << endl;
                        std::cout << "m_sh_csi_weight[" << i << "] of this PE-step:" << m_sh_csi_weight[i] << std::endl;
                        #endif

                        SEQ_TRACE_WARN.print("HEHE4");
                        for (int m = 0; m < 9; m++)
                        {
                            m_adc1[m].getMDH().setClin(static_cast<uint16_t>(i));
                        }

                        sNewRotationMatrix = m_fov->getROT_MATRIX();
                        dRotAngle          = ((2.0 * M_PI) / Nsh) * i;
                        ZRotMatrix[0][0]   = cos(dRotAngle);
                        ZRotMatrix[0][1]   = -sin(dRotAngle);
                        ZRotMatrix[0][2]   = 0;
                        ZRotMatrix[1][0]   = sin(dRotAngle);
                        ZRotMatrix[1][1]   = cos(dRotAngle);
                        ZRotMatrix[1][2]   = 0;
                        ZRotMatrix[2][0]   = 0;
                        ZRotMatrix[2][1]   = 0;
                        ZRotMatrix[2][2]   = 1;
                        sExistingMatrix    = m_fov->getROT_MATRIX();
                        // Matrix Multiplication
                        for (lJ = 0; lJ < 3; lJ++)
                        {
                            for (lK = 0; lK < 3; lK++)
                            {
                                sNewRotationMatrix.dMat[lJ][lK] = 0.0;
                                for (lL = 0; lL < 3; lL++)
                                {
                                    sNewRotationMatrix.dMat[lJ][lK]
                                        += sExistingMatrix.dMat[lJ][lL] * ZRotMatrix[lL][lK];
                                }
                            }
                        }

                        SEQ_TRACE_WARN.print("HEHE5");

                        //if (rMrProt.averages() > j) // measure if this additional scan is really necessary
                        //{
                        //    if (!(m_encod_sl.prepAmplitude(m_asymAmp)) || !(m_encod_sl.check()))
                        //        return m_encod_sl.getNLSStatus();

                        //    // if (!(m_encod_ro.prepAmplitude(
                        //    //         m_d_1st_csi_grad_offset + m_sh_1st_csi_addr[i] * m_d_1st_csi_grad_step))
                        //    //     || !(m_encod_ro.check()))
                        //    //     return m_encod_ro.getNLSStatus();

                        //    // if (!(m_encod_ph.prepAmplitude(
                        //    //         m_d_2nd_csi_grad_offset + m_sh_2nd_csi_addr[i] * m_d_2nd_csi_grad_step))
                        //    //     || !(m_encod_ph.check()))
                        //    //     return m_encod_ph.getNLSStatus();

                            for (int m = 0; m < 9; m++)
                            {
                                m_adc1[m].getMDH().setCacq(static_cast<uint16_t>(j)); // averages
                                // m_adc1[m].getMDH().setCseg(static_cast<uint16_t>(
                                //     m_sh_1st_csi_addr[i] + d1)); // Cseg must be between 0 and
                                //     kSpace().phaseEncodingLines() - 1
                                // m_adc1[m].getMDH().setClin(static_cast<uint16_t>(
                                //     m_sh_2nd_csi_addr[i] + d2)); // Clin must be between 0 and
                                //     kSpace().phaseEncodingLines() - 1
                                // m_adc1[m].getMDH().setCpar(static_cast<uint16_t>(
                                //     m_sh_3rd_csi_addr[i] + d3)); // Cpar must be between 0 and kSpace().partitions()
                                //     - 1
                                SEQ_TRACE_WARN.print("HEHE6");
                                m_adc1[m].getMDH().setCeco(l); // preparing temporal interleaves CF
                                SEQ_TRACE_WARN.print("HEHE7");
                                // flags for extracting time-stamps
                                m_adc1[m].getMDH().setFirstScanInSlice(!i && !j && !l);
                                SEQ_TRACE_WARN.print("HEHE8");
                                m_adc1[m].getMDH().setLastScanInSlice(
                                    i == (static_cast<size_t>(Nsh) - 1) && j == (rMrProt.averages() - 1) && l == (n_ti - 1));
                                SEQ_TRACE_WARN.print("HEHE9");
                                m_adc1[m].getMDH().setLastScanInConcat(
                                    i == (static_cast<size_t>(Nsh) - 1)
                                    && j == (rMrProt.averages() - 1)); // flag for concatenations
                                m_adc1[m].getMDH().setLastScanInMeas(
                                    i == (static_cast<size_t>(Nsh) - 1)
                                    && j == (rMrProt.averages() - 1)); // end of scan
                            }
                            // realize triggering
                            rMrProt.physiology().getPhysioMode(
                                m_firstSignal, m_firstMethod, m_secondSignal, m_secondMethod);

                            if (m_firstMethod
                                == SEQ::METHOD_TRIGGERING /* && m_eTriggerMode == SEQ::TRIGGER_STANDARD */)
                            {
                                long l1 = rMrProt.physiology().triggerDelay(m_firstSignal);
                                long l2 = getMinSBBECGFillTimeRunDelay(
                                    *m_firstHalt); // retrieve min duration of SBBECGFillTimeRun
                                const long triggerDelay = l1 > l2 ? l1 : l2;
                                const long delayPriorToTrigger
                                    = 0; // no gradient ramps must fall within the trigger event block

                                CheckStatusPR(
                                    lStatus = fSBBECGFillTimeRun(&m_firstHalt, triggerDelay, delayPriorToTrigger),
                                    "fSBBECGFillTimeRun");

                                // for estimation  of residual measurement time
                                lCurrKernelCalls++;
                                if (!(lCurrKernelCalls % m_lKernelCallsPerRelevantSignal))
                                {
                                    for (int m = 0; m < 9; m++)
                                    {
                                        m_adc1[m].setRelevantForMeasTime();
                                    }
                                }
                            }

                            // realize off-centre FoV positions by incrementing the phase of the excitation pulses from
                            // step to step
                            double excit_phase = SpecUtils::calcCSIPhase(
                                m_sh_1st_csi_addr[i],
                                m_d_read_pos,
                                rMrProt.sliceSeries().front().readoutFOV(),
                                m_sh_2nd_csi_addr[i],
                                m_d_phase_pos,
                                rMrProt.sliceSeries().front().phaseFOV(),
                                m_sh_3rd_csi_addr[i],
                                m_d_slice_pos,
                                rMrProt.sliceSeries().front().thickness());
                            // double aqc_phase = 0.0;

                            if ((m_sh_1st_csi_addr[i] ^ m_sh_2nd_csi_addr[i] ^ m_sh_3rd_csi_addr[i] ^ j) & 1)
                            {
                                excit_phase += 180.0;
                                // aqc_phase = 180.0;
                            }
                            // else
                            // aqc_phase = 0.0;

                            // m_ph_s_exc.increasePhase(excit_phase);
                            // m_ph_n_exc.decreasePhase(excit_phase);
                            // m_ph_s_adc.increasePhase(aqc_phase);
                            // m_ph_n_adc.decreasePhase(aqc_phase);
                            SEQ_TRACE_WARN.print("HEHE100");
                            lStatus = runKernel(rMrProt, rSeqLim, rSeqExpo, KERNEL_IMAGE);
                            SEQ_TRACE_WARN.print("HEHE101");
                            CheckStatusPR(lStatus, "runKernel");
                            SEQ_TRACE_WARN.print("HEHE102");

                            // undo phase cycling
                            // m_ph_s_exc.decreasePhase(excit_phase);
                            // m_ph_n_exc.increasePhase(excit_phase);
                            // m_ph_s_adc.decreasePhase(aqc_phase);
                            // m_ph_n_adc.increasePhase(aqc_phase);

                        // } // if ( m_sh_csi_weight[i] >= j)
                    }
                }     // PE steps
            }         // averages
        }             // if (rMrProt.kSpace().averagingMode() == SEQ::OUTER_LOOP)    // acquisition outside lines loop

        if (k < (n_meas - 1))
        {
            CheckStatusPR(lStatus = fSBBMeasRepetDelaysRun(rMrProt, rSeqLim, rSeqExpo, k), "fSBBMeasRepetDelaysRun");
        }
    } // end repetition loop

    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ORIGIN_fSEQRunFinish, 0, 0, 0, 0, 0);

    return lStatus;
}

NLS_STATUS Csi_fid::runKernel(
    MrProt&  rMrProt,
    SeqLim&  rSeqLim,
    SeqExpo& rSeqExpo,
    long     lKernelMode,
    long /*lSlice*/,
    long /*lPartition*/,
    long /*lLine*/)
{
    NLS_STATUS lStatus = MRI_SEQ_SEQU_NORMAL;

    MARKER_SEQ_TRACE_DEBUG(99).print("lKernelMode: %ld", lKernelMode);

    // to pass the MDH to the ICE prg.
    for(int m=0;m<9;m++){
    m_adc1[m].getMDH().addToEvalInfoMask(MDH_ONLINE);
    }

    // open this event block
    long loops_per_block= 85; // CF to be changed later
    long spectral_points= 765;
    int n_blocks=spectral_points/loops_per_block;
    fRTEBInit(sNewRotationMatrix, true);
    // it is annoying that opening of an event block is always connected to
    // calculating the rotation matrix which needs to be calculated only once in
    // single slice sequences

    // this timing schemes requires these pre-conditions to be met:
    // - flat top time of exc.gradient > duration of exc. RF puls
    // - ramp down time spoiler == ramp up time slice selection
    // - flat top duration slice selection gradients == RF pulse duration
    // - simultaneously applied spoiling gradients need to be of equal duration

    /************************************* S E Q U E N C E   T I M I N G *************************************/
    /*            Start Time       |    NCO    |  SRF  |  ADC  |    Gradient Events    | Sync                */
    /*              (usec)         |   Event   | Event | Event | phase | read  | slice | Event               */
    /*fRTEI(                       ,           ,       ,       ,       ,       ,       ,        );   [ Clock]*/
    /*********************************************************************************************************/

    fRTEI(0, 0, 0, 0, 0, 0, 0, &m_osc1);
    fRTEI(SysProperties::getCoilCtrlLead(), 0, 0, 0, 0, 0, 0, 0);

    mSEQTest(
        rMrProt,
        rSeqLim,
        rSeqExpo,
        (lKernelMode == KERNEL_CHECK) ? RTEB_ORIGIN_fSEQCheck : RTEB_ORIGIN_fSEQRunKernel,
        10,
        0,
        0,
        0,
        0);

    lStatus = fRTEBFinish();
    if (!MrSucceeded(lStatus))
        return lStatus;

    // Execute Building Blocks if selected:
    //
    //   Water Suppression
    //
    if (!m_sWaterSat.run(rMrProt, rSeqLim, rSeqExpo, nullptr))
        return m_sWaterSat.getNLSStatus();

    fRTEBInit(sNewRotationMatrix, true);

    // NOE
    long lT = 0;
    if (rMrProt.getsSpecPara().getlNOEType() == SEQ::NOE_RECTANGULAR)
    {
        for (int32_t i = 0; i < rMrProt.getsSpecPara().getlNOEs(); i++)
        {
            fRTEI(lT, &m_ph_s_noe, &m_rf_rect_noe, /*A*/ 0, 0, 0, 0, 0);
            lT += m_rf_rect_noe.getDuration();
            fRTEI(lT, &m_ph_n_noe, 0, /*A*/ 0, 0, 0, 0, 0);
            lT += rMrProt.getsSpecPara().getlNOEDelay();
        }
    }

    // excitation
    fRTEI(lT, 0, 0, /*A*/ 0, 0, 0, &m_grad_exc, 0);
    fRTEI(lT += (m_grad_exc.getDuration() - m_rf_exc.getDuration()), &m_ph_s_exc, &m_rf_exc, 0, /*A*/ 0, 0, 0, 0);

    // first spoiler with phase encoding balanced
    fRTEI(lT += m_rf_exc.getDuration(), &m_ph_n_exc, 0, /*A*/ 0, 0, 0, 0, 0);
    fRTEI(lT+= m_grad_exc.getRampDownTime(), 0,0,/*A*/0,&m_RosGy_Prephase, &m_RosGx_Prephase,&m_encod_sl,0);
   	//lT+=m_grad_exc.getRampDownTime() + delay_1;
    // fRTEI(lT+=delay_1, 0,0,/*A*/0,&m_sp1_ph,&m_sp1_ro,&m_sp1_sl,0);
    // lT+=m_sp1_sl.getDuration() + m_sp1_sl.getRampDownTime();
    lT += (m_encod_sl.getTotalTime()+ m_lTEFill_us);


	// refocusing
    // fRTEI(lT+=m_sp1_sl.getDuration(), 0,0,/*A*/ 0,0,0,&m_grad_ref,0);
    // fRTEI(lT+=m_grad_ref.getRampUpTime(), &m_ph_s_ref, &m_rf_ref,0,/*A*/0,0,0,0);

    // second spoiler with rosette prephasing encoding balanced

    // fRTEI(lT+= m_rf_ref.getDuration(), &m_ph_n_ref,0,/*A*/0,&m_sp2_ph,&m_sp2_ro,&m_sp2_sl,0);
	//fRTEI(lT+= m_grad_ref.getRampDownTime(), 0,0,/*A*/0,&m_sp2_ph, &m_sp2_ro,&m_sp2_sl,0);
	//lT+=m_sp2_sl.getDuration() + m_sp2_sl.getRampDownTime() + delay_2;
	// lT+=m_sp2_sl.getDuration() + delay_2;

    // DC
    if (rMrProt.getsSpecPara().getlDecouplingType() == SEQ::DECOUPLING_CW)
    { // acquisition with DC

        const long localTime = lT;
        fRTEI(lT - GRAD_RASTER_TIME * 2, 0, 0, 0, 0, 0, 0, &m_txrfon);
        // Osc bit prevents unnecessary switch events in the "Tx cassette"

        for (long i = 0; i < m_ln_dc_pulses; i++)
        {
            fRTEI(localTime + i * 50000, &m_ph_s_dc, &m_rf_dc, /*A*/ 0, 0, 0, 0, 0);
            fRTEI(localTime + (i + 1) * 50000 - 500, &m_ph_n_dc, 0, /*A*/ 0, 0, 0, 0, 0);
        }
    }

    if (rMrProt.getsSpecPara().getlDecouplingType() == SEQ::DECOUPLING_WALTZ4)
    { // WALTZ 4

        long localTime = lT;

        fRTEI(lT - GRAD_RASTER_TIME * 2, 0, 0, 0, 0, 0, 0, &m_txrfon);
        // Osc bit prevents unnecessary switch events in the "Tx casette"

        for (long i = 0; i < m_ln_dc_pulses; i++)
        {
            fRTEI(localTime, &m_ph_s_dc, &m_rf_waltz, 0, 0, 0, 0, 0);
            localTime += m_rf_waltz.getDuration();
            fRTEI(localTime, &m_ph_s_dc, 0, 0, 0, 0, 0, 0); // since the phase is always nullptr, the negate event
                                                            // m_ph_n_dc is NOT needed
            localTime += m_lPulsePauseRequired;             // delay between pulses must be longer than 100us + 50us
        }
    }

    // acquisition, immediately after the last gradient was switched off
    // fRTEI(m_adc1.getStartTime(), &m_ph_s_adc, 0, &m_adc1, 0, 0, 0, 0);
    // fRTEI(m_adc1.getStartTime() + m_adc1.getRoundedDuration(), &m_ph_n_adc, 0, 0, 0, 0, 0, 0);

    // Rosette ramp up gradients
	
    if (m_ISISalternator == 0)
    {
        fRTEI(lT, 0, 0, 0, &m_RosGy_rampup, &m_RosGx_rampup, 0, 0);
        cout << "Did not alternate." << endl;
    }
    else
    {
        SEQ_TRACE_WARN.print("dalay_3: %ld", delay_3);
        SEQ_TRACE_WARN.print("HOHOHO");
        fRTEI(lT += dummy_delay_3, 0, 0, 0, &m_RosGy_rampup, &m_RosGx_rampup, 0, 0);
        cout << "Alternated!" << endl;
    }
    
    for(int m=0;m<9;m++)
    { 
		m_adc1[m].getMDH().setCpar(0);
    }
    SEQ_TRACE_WARN.print("m_ISISalternator: %ld", m_ISISalternator);
	fRTEI(lT += m_RosGx_rampup.getDuration(), &m_ph_s_adc,            0,   &m_adc1[0],         &m_RosGy,            &m_RosGx, 0 , 0);
    // m_adc1.getMDH().setCseg( 0 );
	//m_adc1.getMDH().setCpar(0);

	fRTEI(lT+=m_RosGy.getDuration(), 0,0,0,0,0,0,0 );
    // mSEQTest(rMrProt, rSeqLim, rSeqExpo, ulTestIdent    , 10, 0/*lLine*/, 0/*lSliceIndex*/, 0, 0) ;

    fRTEBFinish();
	
	for (int i=1; i<n_blocks; i++)
   	{ 
		fRTEBInit( sNewRotationMatrix, true); 
		m_adc1[i].getMDH().setCpar(static_cast<uint16_t>(i));
		fRTEI(0, &m_ph_s_adc,            0,   &m_adc1[i],         &m_RosGy,            &m_RosGx, 0 , 0);
		fRTEI(m_RosGy.getDuration(), 0,0,0,0,0,0,0 );
		lT +=m_RosGy.getDuration();
        // mSEQTest(rMrProt, rSeqLim, rSeqExpo, ulTestIdent    , 10, 0/*lLine*/, 0/*lSliceIndex*/, 0, 0) ;
		fRTEBFinish();
	}

	fRTEBInit( sNewRotationMatrix, true);


    // final spoiling
    // fRTEI(
    //     fSDSRoundUpGRT(m_adc1.getStartTime() + m_adc1.getRoundedDuration()),
    //     0,
    //     0,
    //     0,
    //     &m_finsp_ph,
    //     &m_finsp_ro,
    //     &m_finsp_sl,
    //     0);

    fRTEI(0, &m_ph_n_adc,0,0, /*&m_finsp_ph*/ &m_RosGy_rampdown, /*&m_finsp_ro*/  &m_RosGx_rampdown, &m_finsp_sl, 0 );
    long last_block_dur=(m_finsp_sl.getDuration() + m_finsp_sl.getRampDownTime()+m_RosGy_rampdown.getDuration());
	lT+=last_block_dur;

    // TR fill
    // fRTEI(
    //     fSDSRoundUpGRT(m_adc1.getStartTime() + m_adc1.getRoundedDuration()) + m_finsp_ph.getTotalTime() + m_lTRFill,
    //     0,
    //     0,
    //     0,
    //     0,
    //     0,
    //     0,
    //     0);
    fRTEI(last_block_dur+(rMrProt.tr()[0] - lT-m_sWaterSat.getDurationPerRequest()-SysProperties::getCoilCtrlLead()), 0,0,0,0,0,0,0);
    SEQ_TRACE_WARN.print("HEHEHE1");

    // do testing and close the event block
    mSEQTest(rMrProt, rSeqLim, rSeqExpo, RTEB_ClockCheck, 10, 0 /*lLine*/, 0 /*lSliceIndex*/, 0, 0);
    mSEQTest(
        rMrProt,
        rSeqLim,
        rSeqExpo,
        (lKernelMode == KERNEL_CHECK) ? RTEB_ORIGIN_fSEQCheck : RTEB_ORIGIN_fSEQRunKernel,
        10,
        0 /*lLine*/,
        0 /*lSliceIndex*/,
        0,
        0);

    SEQ_TRACE_WARN.print("HEHEHE2");
    SEQ_TRACE_WARN.print("dalay_3: %ld", delay_3);
    CheckStatusPR(lStatus = fRTEBFinish(), "fRTEBFinish [*0010*]");
    SEQ_TRACE_WARN.print("HEHEHE3");
    return lStatus;
}

void Csi_fid::createUI()
{
#ifdef WIN32
    delete m_pUI;
    m_pUI = new Csi_fidUI();
#endif
}

#ifdef WIN32
SpecUI* Csi_fid::getUI()
{
    return m_pUI;
}
#endif

long Csi_fid::getTRNeededSpectro()
{
    return m_lTRNeededSpectro;
}

void Csi_fid::setTRNeededSpectro(long value)
{
    m_lTRNeededSpectro = value;
}
