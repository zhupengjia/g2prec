// -*- C++ -*-

/* class G2PRec
 * The real class to do the reconstruction.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Feb 2014, C. Gu, Modified for G2PRec.
//

#ifndef G2P_REC_H
#define G2P_REC_H

#include "G2PAppBase.hh"

class G2PDrift;

class G2PRec : public G2PAppBase {
public:
    G2PRec();
    virtual ~G2PRec();

    virtual int Process(const float* V5bpm_bpm, const float* V5bpmave_bpm, const double* V5tp_tr, double* V5corr_tr, double* V5rec_tr, double* V5rec_lab);

    void SetBeamEnergy(double e);
    void SetHRSMomentum(double p);

protected:
    virtual int Initialize();
    virtual void Clear(Option_t* /*option*/ = "");

    void TransBPM2Tr(const double* V5_bpm, double* V5_tr);

    void GetEffBPM(const double* V2_tr, double* V2eff_tr);
    void ExtTgtCorr(double xbpm, double ybpm, const double* V5tp_tr, double* V5tpcorr_tr);

    virtual int Configure();

    double fBeamEnergy;
    double fHRSAngle;
    double fHRSMomentum;
    double fFieldRatio;

    double fDriftLimit;

    double fSieveZ;
    double fRecZ;

    double fExtTgtCorrT;
    double fExtTgtCorrP;
    double fExtTgtCorrD;

    double fFitPars[2][3];

    double fV5bpm_bpm[5];
    double fV5bpm_tr[5];

    double fV5bpmave_bpm[5];
    double fV5bpmave_tr[5];

    double fV5tpmat_tr[5];
    double fV5sieveproj_tr[5];

    double fV5rec_tr[5];
    double fV5rec_lab[5];

    G2PDrift* pDrift;

private:
    ClassDef(G2PRec, 1)
};

#endif