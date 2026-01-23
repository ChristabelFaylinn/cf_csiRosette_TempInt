//-----------------------------------------------------------------------------
// <copyright file="csi_fid.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
// <description>CSI FID sequence, no in-plane VoI selection</description>
//-----------------------------------------------------------------------------

#pragma once

#include "MrGlobalDefinitions/MrResult.h"
#include "MrImagingFW/libSBBFW/StdSeqIF.h"
#include "MrImagingFW/libSBBFW/libSBBFW.h"
#include "MrMeasSrv/SeqIF/libRT/sFREQ_PHASE.h"
#include "MrMeasSrv/SeqIF/libRT/sGRAD_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sREADOUT.h"
#include "MrMeasSrv/SeqIF/libRT/sRF_PULSE.h"
#include "MrMeasSrv/SeqIF/libRT/sSLICE_POS.h"
#include "MrMeasSrv/SeqIF/libRT/sSYNC.h"
#include "MrProtSrv/Domain/CoreNative/MrSpectroscopy.h"
#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrProt.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/SeqIF/SeqExpo.h"
#include "MrSpecAcq/SpecSeq/SpecSeq.h"
#include "MrSpecAcq/libSpecSBB/SBBWatSup.h"

#include <array>
#include <vector>

#ifdef BUILD_SEQU
#define __OWNER
#endif
#include "MrGlobalDefinitions/ImpExpCtrl.h"

class MrProt;
class SeqLim;
class SeqExpo;

namespace SEQ_NAMESPACE
{
class Csi_fidUI;

class __IMP_EXP Csi_fid : public SpecSeq
{
  public:
    Csi_fid() = default;
    virtual ~Csi_fid();

    Csi_fid(const Csi_fid& right) = delete;
    Csi_fid& operator=(const Csi_fid& right) = delete;
    Csi_fid(Csi_fid&& right)                 = delete;
    Csi_fid& operator=(Csi_fid&& right) = delete;

    NLSStatus  initialize(SeqLim& rSeqLim) override;
    NLSStatus  prepare(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;
    NLSStatus  check(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo, SEQCheckMode* pSEQCheckMode) override;
    NLSStatus  run(MrProt& rMrProt, SeqLim& rSeqLim, SeqExpo& rSeqExpo) override;
    NLS_STATUS runKernel(
        MrProt&  rMrProt,
        SeqLim&  rSeqLim,
        SeqExpo& rSeqExpo,
        long     lKernelMode,
        long     lSlice     = 0,
        long     lPartition = 0,
        long     lLine      = 0) override;

#ifdef WIN32
    virtual SpecUI* getUI();
#endif

    long getTRNeededSpectro() override;
    void setTRNeededSpectro(long value);

  protected:
    static const int32_t kWaltzMaxNoSamples{4096};

    // Slice position information (rotation matrix and voxel position)
    sSLICE_POS m_fov;
    sSLICE_POS m_voi;

    // Sync Bits
    sSYNC m_osc1{"osc1"};
    sSYNC m_txrfon{"txrfon"};

    // RF-Pulses and NO-Events
    sRF_PULSE_EXT m_rf_exc{"rf_exc"};

    sFREQ_PHASE m_ph_s_exc{"ph_s_exc"};
    sFREQ_PHASE m_ph_n_exc{"ph_n_exc"};

    // added -------------------------
    sRF_PULSE_EXT m_rf_ref{"rf_ref"};
    sFREQ_PHASE m_ph_s_ref{"ph_s_ref"};
    sFREQ_PHASE m_ph_n_ref{"ph_n_ref"};

    // NOE pulse
    sRF_PULSE_RECT m_rf_rect_noe{"rf_rect_noe"};

    sFREQ_PHASE m_ph_s_noe{"ph_s_noe"};
    sFREQ_PHASE m_ph_n_noe{"ph_n_noe"};

    // DC pulse
    sRF_PULSE_RECT m_rf_dc{"rf_dc"}; // used also for the 180 deg. WALTZ pulse

    sFREQ_PHASE m_ph_s_dc{"ph_s_dc"};
    sFREQ_PHASE m_ph_n_dc{"ph_n_dc"};

    sRF_PULSE_ARB m_rf_waltz{"rf_waltz"};

    std::array<sSample, kWaltzMaxNoSamples> m_samples{};

    // measurement time
    long m_lKernelCallsPerRelevantSignal{1};

    // Physio
    sSYNC_PHYSIO1_HALT m_firstHalt;
    SEQ::PhysioSignal  m_firstSignal{SEQ::SIGNAL_NONE};
    SEQ::PhysioMethod  m_firstMethod{SEQ::METHOD_NONE};
    SEQ::PhysioSignal  m_secondSignal{SEQ::SIGNAL_NONE};
    SEQ::PhysioMethod  m_secondMethod{SEQ::METHOD_NONE};

    // Readout Events
    // sREADOUT m_adc1{"adc1"};

    sREADOUT m_adc1[12];
    
    sFREQ_PHASE m_ph_s_adc{"m_ph_s_adc"};
    sFREQ_PHASE m_ph_n_adc{"m_ph_n_adc"};

    /////[////////////////////
    // Rosette Gradients
    ////////////////////////

    sGRAD_PULSE_ARB m_RosGx{"m_Ros_Gx"};
    sGRAD_PULSE_ARB m_RosGy{"m_Ros_Gy"};
    sGRAD_PULSE_ARB m_RosGx_rampup{"m_Ros_Gx_rampup"};
    sGRAD_PULSE_ARB m_RosGy_rampup{"m_Ros_Gy_rampup"};
    sGRAD_PULSE_ARB m_RosGx_rampdown{"m_Ros_Gx_rampdown"};
    sGRAD_PULSE_ARB m_RosGy_rampdown{"m_Ros_Gy_rampdown"};
    sGRAD_PULSE     m_RosRewx{"m_RosRewx"};
    sGRAD_PULSE     m_RosRewy{"m_RosRewy"};   
    sGRAD_PULSE     m_RosGx_Prephase{"m_RosGx_Prephase"};
    sGRAD_PULSE     m_RosGy_Prephase{"m_RosGy_Prephase"};

    // Gradient Pulses
    // slice selection pulses
    sGRAD_PULSE m_grad_exc{"grad_exc"}; // gradient during excitation
    sGRAD_PULSE m_grad_ref{"grad_ref"};
    
    // phase encoding pulses
    sGRAD_PULSE m_encod_sl{"encod_sl"};
    sGRAD_PULSE m_encod_ph{"encod_ph"};
    sGRAD_PULSE m_encod_ro{"encod_ro"};

    // ------------------------------------------------------------------------------
    /// TE fill time in us
    // ------------------------------------------------------------------------------
    long m_lTEFill_us{0};
    long delay_1{0};
    long delay_2{0};
    long delay_3{0}; // CF for temp int: delay between temp interleaf
    long dummy_delay_3{0};

    /////////////////
    // spoiler pulses
    /////////////////

    sGRAD_PULSE     m_sp1_ph{"sp1_ph"}; // spoiler 1 - phase
    sGRAD_PULSE     m_sp1_ro{"sp1_ro"}; // spoiler 1 - readout
    sGRAD_PULSE     m_sp1_sl{"sp1_s1"}; // spoiler 1 - slice

    sGRAD_PULSE     m_sp2_ph{"sp2_ph"}; // spoiler 2 - phase
    sGRAD_PULSE     m_sp2_ro{"sp2_ro"}; // spoiler 2 - readout
    sGRAD_PULSE     m_sp2_sl{"sp2_s1"}; // spoiler 2 - slice

    // spoiler pulses
    sGRAD_PULSE m_finsp_sl{"finsp_sl"}; // final spoiler - slice
    sGRAD_PULSE m_finsp_ro{"finsp_ro"}; // final spoiler - readout
    sGRAD_PULSE m_finsp_ph{"finsp_ph"}; // final spoiler - phase

    // MR Spectroscopy Sequence Functionality
    long m_lSamplesBeforeEcho{0};

    // CSI variables
    std::vector<int32_t> m_sh_1st_csi_addr;
    std::vector<int32_t> m_sh_2nd_csi_addr;
    std::vector<int32_t> m_sh_3rd_csi_addr;
    std::vector<int32_t> m_sh_csi_weight;


    double m_d_1st_csi_grad_step{0.0};
    double m_d_2nd_csi_grad_step{0.0};
    double m_d_3rd_csi_grad_step{0.0};
    double m_d_1st_csi_grad_offset{0.0};
    double m_d_2nd_csi_grad_offset{0.0};
    double m_d_3rd_csi_grad_offset{0.0};
    double m_d_read_pos{0.0};
    double m_d_phase_pos{0.0};
    double m_d_slice_pos{0.0};

    double m_asymAmp{0.0};  // amplitude defined for implementing asymmmetric RF pulse

    // DC variables
    long m_ln_dc_pulses{0};
    long m_lPulsePauseRequired{0};

    // Timing variables
    long m_lTRNeededSpectro{0};
    long m_lTRFill{0};
    long m_ISISalternator{1};

    // SBBs
    SeqBuildBlockWatSup m_sWaterSat{nullptr};

#ifdef WIN32
    Csi_fidUI* m_pUI{nullptr};
#endif

    virtual void createUI();
};
} // namespace SEQ_NAMESPACE
