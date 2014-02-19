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

    virtual int Process(const float* V5bpm_bpm, const double* V5tp_tr, double* V5rec_tr, double* V5rec_lab);

    void SetHRSMomentum(double p);

protected:
    virtual int Initialize();
    virtual void Clear();

    double GetEffBPMX();

    virtual int Configure();

    double fHRSAngle;
    double fHRSMomentum;
    double fFieldRatio;

    double fSieveZ;
    double fRecZ;

    double fV5bpm_bpm[5];

    double fV5tpmat_tr[5];
    double fV5sieveproj_tr[5];

    double fV5tprec_tr[5];
    double fV5tprec_lab[5];

    G2PDrift* pDrift;

private:
    ClassDef(G2PRec, 1)
};

#endif