//-----------------------------------------------------------------------------
// <copyright file="csi_fid_UI.h" company="Siemens Healthcare GmbH">
//   Copyright (C) Siemens Healthcare GmbH, 2020. All Rights Reserved. Confidential.
// </copyright>
//-----------------------------------------------------------------------------

#pragma once

#ifndef csi_fid_UI_h
#define csi_fid_UI_h

#ifdef WIN32

#include "MrProtSrv/Domain/CoreNative/SeqLim.h"
#include "MrProtSrv/Domain/MrProtocol/libUICtrl/UICtrl.h"
#include "MrSpecAcq/SpecSeq/SpecUI.h"

class SeqLim;

namespace SEQ_NAMESPACE
{
class Csi_fidUI : public SpecUI
{
  public:
    Csi_fidUI()          = default;
    virtual ~Csi_fidUI() = default;

    Csi_fidUI(const Csi_fidUI& right) = delete;
    Csi_fidUI& operator=(const Csi_fidUI& right) = delete;
    Csi_fidUI(Csi_fidUI&& right)                 = delete;
    Csi_fidUI& operator=(Csi_fidUI&& right) = delete;

    //  --------------------------------------------------------------------------
    //
    //  Name        : Csi_fidUI::registerUI
    //
    //  Description :
    /// \brief        This function initializes the UI functions and
    ///                registers all given set / get / Solve - handlers
    ///
    ///               It can be executed on the measurement system, too, but is empty there.
    ///
    ///               On the host, it executes these steps
    ///               - Declaration of pointers to UI classes
    ///               - Registration of overloaded set value handlers
    ///
    ///               It returns an NLS status
    ///
    NLS_STATUS addFunctionalityToRegisterUI(SeqLim& rSeqLim) override;

    UI_ELEMENT_DOUBLE m_CsiVoiFovRead;
    UI_ELEMENT_DOUBLE m_CsiVoiFovPhase;
    UI_ELEMENT_DOUBLE m_CsiReadOutFov;
    UI_ELEMENT_DOUBLE m_CsiPhaseFov;
    UI_ELEMENT_DOUBLE m_DecouTotalDur;
};
} // namespace SEQ_NAMESPACE

#endif // WIN32
#endif // csi_fid_UI_h
