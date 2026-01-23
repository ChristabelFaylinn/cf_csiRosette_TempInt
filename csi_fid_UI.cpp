//-----------------------------------------------------------------------------
// <copyright file="csi_fid_UI.cpp" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
//-----------------------------------------------------------------------------

#ifdef WIN32

#include "MrSpecAcq/cf_csiRosette_TempInt/csi_fid_UI.h"

#include "MrMeasSrv/SeqIF/Sequence/sequmsg.h"
#include "MrProtSrv/Domain/MrProtData/MrProt/MrSliceGroup.h"
#include "MrProtSrv/Domain/MrProtocol/UILink/MrStdNameTags.h"
#include "MrSpecAcq/cf_csiRosette_TempInt/csi_fid.h"

#ifndef SEQ_NAMESPACE
#error SEQ_NAMESPACE not defined
#endif
using namespace SEQ_NAMESPACE;

namespace Csi_fidUINS
{
//  ----------------------------------------------------------------------
//
//  Name        :  getSeq
//
//  Description :
/// \brief         Returns the pointer to the sequence Csi_fid
//
//  Return      :  Csi_fid*
//
//  ----------------------------------------------------------------------
Csi_fid* getSeq(MrUILinkBase* const pThis)
{
    return static_cast<Csi_fid*>(pThis->sequence().getSeq());
}

// * -------------------------------------------------------------------------- *
// * -------------------------------------------------------------------------- *
// *                                                                            *
// *                    Definition of get Limits handler                        *
// *                                                                            *
// * -------------------------------------------------------------------------- *
// * -------------------------------------------------------------------------- *

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// set number of measurements limit to Hard Limit (currently 8) for 2D and to 4 for 3D    //overload
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

static bool fmyUILinkMeasurementsGetLimits(
    LINK_LONG_TYPE* const pThis, std::vector<MrLimitLong>& rLimitVector, uint32_t& rulVerify, int32_t /*pos*/)
{
    rLimitVector.resize(1);
    rulVerify                         = LINK_LONG_TYPE::VERIFY_BINARY_SEARCH;
    const ParLim<int32_t>& _seqLimits = pThis->seqLimits().getRepetitions();

    const long lMax = 4;

    if (pThis->prot().getsKSpace().getucDimension() == SEQ::DIM_3)
    {
        return rLimitVector[0].setEqualSpaced((1 + _seqLimits.getMin()), lMax, _seqLimits.getInc());
    }

    // case 2D : hard limits
    return rLimitVector[0].setEqualSpaced((1 + _seqLimits.getMin()), (1 + _seqLimits.getMax()), _seqLimits.getInc());
}

// * -------------------------------------------------------------------------- *
// *                                                                            *
// *            Limit Handlers                                             *
// *                                                                            *
// * -------------------------------------------------------------------------- *

void Init_GetLimitHandler_NofMeasurementLimit(SeqLim& rSeqLim)
{
    // overload
    // register GetLimitHandler Number of Measurements
    LINK_LONG_TYPE* pNofMeasurementLimit = _search<LINK_LONG_TYPE>(rSeqLim, MR_TAG_MEASUREMENTS);
    if (pNofMeasurementLimit)
    {
        pNofMeasurementLimit->registerGetLimitsHandler(fmyUILinkMeasurementsGetLimits);
    }
}

//---------------------------------------------------------------------------
// for TD/ pulse duration:
//---------------------------------------------------------------------------

#ifdef SUPPORT_MULTISPEC // only for xFID

void Init_AvailableHandler_ExcPulseDur(SeqLim& rSeqLim)
{
    // prototype:
    bool fSAvailableExcPulseDuration(LINK_DOUBLE_TYPE* const pThis, long lIndex);

    LINK_DOUBLE_TYPE::PFctIsAvailable pExcPulseDurationIsAvailableHandler;
    // Declaration of pointers to UI classes
    LINK_DOUBLE_TYPE* pExcPulseDuration = _searchElm<LINK_DOUBLE_TYPE>(rSeqLim, MR_TAG_TD); // MR_TAG_TE

    // Registration of is available handler
    if (pExcPulseDuration)
    {
        pExcPulseDurationIsAvailableHandler
            = pExcPulseDuration->registerIsAvailableHandler(fSAvailableExcPulseDuration);
    }
}

#endif

unsigned myDecouplingTotalDurationGetToolTip(LINK_DOUBLE_TYPE* const pThis, char** arg_list, int32_t)
{
    static char tToolTip[12];
    sprintf(tToolTip, "%u ms", pThis->prot().getaulServicePara()[1]);
    arg_list[0] = tToolTip;
    return MRI_STD_STRING; // for hard coded text
}

/*[ Function ****************************************************************\
*
* Name        : myCsiVoiFovReadGetLimits, myCsiVoiFovPhaseGetLimits
*
* Description : VOI read direction will be dimmed, no in-plane VOI in CSI_FID
*
* Return      : bool
*
\****************************************************************************/

static bool myCsiVoiFovReadGetLimits(
    LINK_DOUBLE_TYPE* const /*pThis*/,
    std::vector<MrLimitDouble>& /*rLimitVector*/,
    uint32_t& /*rulVerify*/,
    int32_t /*pos*/)
{
    return false;
}

static bool myCsiVoiFovPhaseGetLimits(
    LINK_DOUBLE_TYPE* const /*pThis*/,
    std::vector<MrLimitDouble>& /*rLimitVector*/,
    uint32_t& /*rulVerify*/,
    int32_t /*pos*/)
{
    return false;
}

static double myCsiFovReadSetValue(LINK_DOUBLE_TYPE* const pThis, double newVal, int32_t /*pos*/)
{
    MrProt rMrProt(pThis->prot());
    if (rMrProt.sliceSeries().empty())
        return 0;

    SliceGroupArray& rGroupArray = rMrProt.sliceGroupList();
    for (long cntr = 0; cntr < rGroupArray.size(); ++cntr)
    {
        rGroupArray[cntr].readoutFOV(newVal);
        rMrProt.getsSpecPara().getsVoI().setdReadoutFOV(newVal);
    }
    return rMrProt.sliceSeries().aFront().readoutFOV();
}

static double myCsiFovPhaseSetValue(LINK_DOUBLE_TYPE* const pThis, double newVal, int32_t /*pos*/)
{
    MrProt rMrProt(pThis->prot());
    if (rMrProt.sliceSeries().empty())
        return 0;

    SliceGroupArray& rGroupArray = rMrProt.sliceGroupList();
    for (long cntr = 0; cntr < rGroupArray.size(); ++cntr)
    {
        rGroupArray[cntr].phaseFOV(newVal);
        rMrProt.getsSpecPara().getsVoI().setdPhaseFOV(newVal);
    }
    return rMrProt.sliceSeries().aFront().phaseFOV();
}

} // namespace Csi_fidUINS

//  --------------------------------------------------------------------------
//
//  Name        : registerUI
//
//  Description :
/// \brief        This method registers all given set / get / Solve - handlers
///
///               It can be executed on the measurement system, too, but is empty there.
///
///
///               It returns an NLS status
///
//  Return      : long
//
//  --------------------------------------------------------------------------

NLS_STATUS Csi_fidUI::addFunctionalityToRegisterUI(SeqLim& rSeqLim)
{
#ifdef WIN32

    // register the redefined functions specified above with the UI
    m_CsiVoiFovRead.registerGetLimitsHandler(rSeqLim, MR_TAG_CSI_VOI_FOV_READ, Csi_fidUINS::myCsiVoiFovReadGetLimits);
    m_CsiVoiFovPhase.registerGetLimitsHandler(
        rSeqLim, MR_TAG_CSI_VOI_FOV_PHASE, Csi_fidUINS::myCsiVoiFovPhaseGetLimits);
    m_CsiReadOutFov.registerSetValueHandler(rSeqLim, MR_TAG_CSI_READOUT_FOV, Csi_fidUINS::myCsiFovReadSetValue);
    m_CsiPhaseFov.registerSetValueHandler(rSeqLim, MR_TAG_CSI_PHASE_FOV, Csi_fidUINS::myCsiFovPhaseSetValue);
    m_DecouTotalDur.registerToolTipHandler(
        rSeqLim, MR_TAG_DECOUPLING_TOTAL_DURATION, Csi_fidUINS::myDecouplingTotalDurationGetToolTip);

    // UI-utilities defined in spectro_ui.cpp

    Csi_fidUINS::Init_GetLimitHandler_NofMeasurementLimit(rSeqLim); // overload GetLimitHandler Number of Measurement

#endif

    return MRI_SEQ_SEQU_NORMAL;
}

#endif // WIN32
