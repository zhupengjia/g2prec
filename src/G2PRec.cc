// -*- C++ -*-

/* class G2PRec
 * The real class to do the reconstruction.
 */

// History:
//   Apr 2013, C. Gu, First public version.
//   Feb 2014, C. Gu, Modified for G2PRec.
//

#include <cstdlib>
#include <cstdio>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TObject.h"

#include "libconfig.h++"

#include "G2PAppBase.hh"
#include "G2PConf.hh"
#include "G2PDrift.hh"

#include "G2PRec.hh"

using namespace std;
using namespace libconfig;

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PRec::G2PRec() :
fBeamEnergy(2.254), fHRSAngle(5.767 * kDEG), fHRSMomentum(2.254), fFieldRatio(0.0), fSieveZ(800.0), fRecZ(0.0), fExtTgtCorrT(0.0), fExtTgtCorrP(0.0), fExtTgtCorrD(0.0)
{
    // Constructor

    memset(fFitPars, 0, sizeof (fFitPars));

    pDrift = new G2PDrift();

    Configure();
    Initialize();
    Clear();
}

G2PRec::~G2PRec()
{
    // Destructor

    delete pDrift;
    pDrift = NULL;
}

int G2PRec::Process(const float* V5bpm_bpm, const float* V5bpmave_bpm, const double* V5tp_tr, double* V5corr_tr, double* V5rec_tr, double* V5rec_lab)
{
    static const char* const here = "Process()";

    fV5bpm_bpm[0] = V5bpm_bpm[0] / 1000.0;
    fV5bpm_bpm[1] = V5bpm_bpm[1];
    fV5bpm_bpm[2] = V5bpm_bpm[2] / 1000.0;
    fV5bpm_bpm[3] = V5bpm_bpm[3];
    fV5bpm_bpm[4] = V5bpm_bpm[4] / 1000.0;

    fV5bpmave_bpm[0] = V5bpmave_bpm[0] / 1000.0;
    fV5bpmave_bpm[1] = 0.0;
    fV5bpmave_bpm[2] = V5bpmave_bpm[2] / 1000.0;
    fV5bpmave_bpm[3] = 0.0;
    fV5bpmave_bpm[4] = fV5bpm_bpm[4];

    fRecZ = fV5bpm_bpm[4];

    TransBPM2Tr(fV5bpm_bpm, fV5bpm_tr);
    TransBPM2Tr(fV5bpmave_bpm, fV5bpmave_tr);

    ExtTgtCorr(fV5bpmave_tr[0], fV5bpmave_tr[2], V5tp_tr, fV5tpmat_tr);

    /*
    double delta_old = 1e38;
    double xeff_tr = 0, yeff_tr = 0;
    while (fabs(fV5tpmat_tr[4] - delta_old) > 1.e-4) {
        delta_old = fV5tpmat_tr[4];
        xeff_tr = GetEffBPM(0);
        yeff_tr = GetEffBPM(1);
        ExtTgtCorr(xeff_tr, yeff_tr, V5tp_tr, fV5tpmat_tr);
    }*/

    double V2[2], V2eff[2];
    V2[0] = fV5bpmave_tr[0];
    V2[1] = fV5bpmave_tr[2];
    GetEffBPM(V2, V2eff);
    ExtTgtCorr(V2eff[0], V2eff[1], V5tp_tr, fV5tpmat_tr);

    V2[0] = fV5bpm_tr[0];
    V2[1] = fV5bpm_tr[2];
    GetEffBPM(V2, V2eff);
    fV5tpmat_tr[0] = V2eff[0];
    fV5tpmat_tr[2] = V2eff[1];

    if (fDebug > 0) {
        Info(here, "bpm_bpm   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_bpm[0], fV5bpm_bpm[1], fV5bpm_bpm[2], fV5bpm_bpm[3], fV5bpm_bpm[4]);
        Info(here, "tpmat_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpmat_tr[0], fV5tpmat_tr[1], fV5tpmat_tr[2], fV5tpmat_tr[3], fV5tpmat_tr[4]);
    }

    for (int i = 0; i < 5; i++) {
        V5corr_tr[i] = fV5tpmat_tr[i];
    }
    V5corr_tr[1] = tan(fV5tpmat_tr[1]);
    V5corr_tr[3] = tan(fV5tpmat_tr[3]);

    Project(fV5tpmat_tr[0], fV5tpmat_tr[2], 0.0, fSieveZ, fV5tpmat_tr[1], fV5tpmat_tr[3], fV5sieveproj_tr[0], fV5sieveproj_tr[2]);
    fV5sieveproj_tr[1] = fV5tpmat_tr[1];
    fV5sieveproj_tr[3] = fV5tpmat_tr[3];
    fV5sieveproj_tr[4] = fV5tpmat_tr[4];

    if (fDebug > 0) {
        Info(here, "sivproj_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieveproj_tr[0], fV5sieveproj_tr[1], fV5sieveproj_tr[2], fV5sieveproj_tr[3], fV5sieveproj_tr[4]);
    }

    if (pDrift->Drift(fV5sieveproj_tr, fSieveZ, fHRSMomentum, fHRSAngle, 0.0, fV5rec_tr) > fDriftLimit) {
        for (int i = 0; i < 5; i++) {
            V5rec_tr[i] = 1.e+38;
            V5rec_lab[i] = 1.e+38;
        }
        return -1;
    }

    TCS2HCS(fV5rec_tr[0], fV5rec_tr[2], 0.0, fHRSAngle, fV5rec_lab[0], fV5rec_lab[2], fV5rec_lab[4]);
    TCS2HCS(fV5rec_tr[1], fV5rec_tr[3], fHRSAngle, fV5rec_lab[1], fV5rec_lab[3]);

    double x[3] = {fV5rec_lab[0], fV5rec_lab[2], fV5rec_lab[4]};
    double p[3] = {fHRSMomentum * (1 + fV5rec_tr[4]) * sin(fV5rec_lab[1]) * cos(fV5rec_lab[3]),
        fHRSMomentum * (1 + fV5rec_tr[4]) * sin(fV5rec_lab[1]) * sin(fV5rec_lab[3]),
        fHRSMomentum * (1 + fV5rec_tr[4]) * cos(fV5rec_lab[1])};

    if (pDrift->Drift(x, p, fRecZ, x, p) > fDriftLimit) {
        for (int i = 0; i < 5; i++) {
            V5rec_tr[i] = 1.e+38;
            V5rec_lab[i] = 1.e+38;
        }
        return -1;
    };

    double temp;
    fV5rec_lab[0] = x[0];
    fV5rec_lab[1] = acos(p[2] / (fHRSMomentum * (1 + fV5rec_tr[4])));
    fV5rec_lab[2] = x[1];
    fV5rec_lab[3] = atan2(p[1], p[0]);
    fV5rec_lab[4] = x[2];
    HCS2TCS(x[0], x[1], x[2], fHRSAngle, fV5rec_tr[0], fV5rec_tr[2], temp);
    HCS2TCS(fV5rec_lab[1], fV5rec_lab[3], fHRSAngle, fV5rec_tr[1], fV5rec_tr[3]);

    if (fDebug > 0) {
        Info(here, "rec_tr    : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rec_tr[0], fV5rec_tr[1], fV5rec_tr[2], fV5rec_tr[3], fV5rec_tr[4]);
        Info(here, "rec_lab   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5rec_lab[0], fV5rec_lab[1], fV5rec_lab[2], fV5rec_lab[3], fV5rec_lab[4]);
    }

    for (int i = 0; i < 5; i++) {
        V5rec_tr[i] = fV5rec_tr[i];
        V5rec_lab[i] = fV5rec_lab[i];
    }
    V5rec_tr[1] = tan(V5rec_tr[1]);
    V5rec_tr[3] = tan(V5rec_tr[3]);

    return 0;
}

void G2PRec::SetBeamEnergy(double e)
{
    static const char* const here = "Configure()";

    fBeamEnergy = e;

    if (fDebug > 0) {
        Info(here, "fBeamEnergy\t= %le", fBeamEnergy);
    }
}

void G2PRec::SetHRSMomentum(double p)
{
    static const char* const here = "Configure()";

    fHRSMomentum = p;

    if (fDebug > 0) {
        Info(here, "fHRSMomentum\t= %le", fHRSMomentum);
    }
}

int G2PRec::Initialize()
{
    //static const char* const here = "Initialize()";

    return 0;
}

void G2PRec::Clear(Option_t* /*option*/)
{
    memset(fV5bpm_bpm, 0, sizeof (fV5bpm_bpm));
    memset(fV5bpm_tr, 0, sizeof (fV5bpm_tr));
    memset(fV5bpmave_bpm, 0, sizeof (fV5bpmave_bpm));
    memset(fV5bpmave_tr, 0, sizeof (fV5bpmave_tr));
    memset(fV5tpmat_tr, 0, sizeof (fV5tpmat_tr));
    memset(fV5sieveproj_tr, 0, sizeof (fV5sieveproj_tr));
    memset(fV5rec_tr, 0, sizeof (fV5rec_tr));
    memset(fV5rec_lab, 0, sizeof (fV5rec_lab));
}

void G2PRec::TransBPM2Tr(const double* V5_bpm, double* V5_tr)
{
    static const char* const here = "TransBPM2Tr()";

    //double x[3] = {V5_bpm[0], V5_bpm[2], V5_bpm[4]};
    double p[3] = {tan(V5_bpm[3]), tan(V5_bpm[1]), 1.0};
    double pp = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    double theta = acos(1.0 / pp);
    double phi = atan2(p[1], p[0]);

    HCS2TCS(V5_bpm[0], V5_bpm[2], V5_bpm[4], fHRSAngle, V5_tr[0], V5_tr[2], V5_tr[4]);
    HCS2TCS(theta, phi, fHRSAngle, V5_tr[1], V5_tr[3]);

    if (fDebug > 2) {
        Info(here, "%10.3e %10.3e %10.3e -> %10.3e %10.3e %10.3e", V5_bpm[0], V5_bpm[2], V5_bpm[4], V5_tr[0], V5_tr[2], V5_tr[4]);
    }
}

void G2PRec::GetEffBPM(const double* V2_tr, double* V2eff_tr)
{
    static const char* const here = "GetEffBPM()";

    V2eff_tr[0] = V2_tr[0];
    V2eff_tr[1] = V2_tr[1];

    if (fFieldRatio > 1e-5) {
        // Fit:
        // (Xbpm_tr-Xeffbpm_tr) vs P
        // ([0]+[1]/x)
        double p = (1 + fV5tpmat_tr[4]) * fHRSMomentum;
        V2eff_tr[0] = V2_tr[0] - (fFitPars[0][0] + (fFitPars[0][1] + fFitPars[0][2] * V2_tr[1]) / p) / 1000;
        V2eff_tr[1] = V2_tr[1] - (fFitPars[1][0] + (fFitPars[1][1] + fFitPars[1][2] * V2_tr[0]) / p) / 1000;
    }

    if (fDebug > 1) {
        Info(here, "effbpm_tr : %10.3e %10.3e", V2eff_tr[0], V2eff_tr[1]);
    }
}

void G2PRec::ExtTgtCorr(double xbpm, double ybpm, const double* V5tp_tr, double* V5tpcorr_tr)
{
    V5tpcorr_tr[0] = V5tp_tr[0];
    V5tpcorr_tr[1] = V5tp_tr[1] + fExtTgtCorrT * xbpm;
    V5tpcorr_tr[2] = V5tp_tr[2];
    V5tpcorr_tr[3] = V5tp_tr[3] + fExtTgtCorrP * ybpm;
    V5tpcorr_tr[4] = V5tp_tr[4] + xbpm / fExtTgtCorrD;
}

int G2PRec::Configure()
{
    static const char* const here = "Configure()";

    if (G2PAppBase::Configure() != 0) return -1;

    if (!gConfig->lookupValue("field.ratio", fFieldRatio))
        Warning(here, "Cannot find setting \"field.ratio\", using default value ......");

    if (!gConfig->lookupValue("hrs.angle", fHRSAngle))
        Warning(here, "Cannot find setting \"hrs.angle\", using default value ......");

    if (!gConfig->lookupValue("sieve.z", fSieveZ))
        Warning(here, "Cannot find setting \"sieve.z\", using default value ......");

    if (!(gConfig->lookupValue("rec.fit.x.p0", fFitPars[0][0])
            && gConfig->lookupValue("rec.fit.x.p1", fFitPars[0][1])
            && gConfig->lookupValue("rec.fit.x.p2", fFitPars[0][2])))
        Warning(here, "Cannot find setting \"rec.beamfit\", using default value ......");

    if (!(gConfig->lookupValue("rec.fit.y.p0", fFitPars[1][0])
            && gConfig->lookupValue("rec.fit.y.p1", fFitPars[1][1])
            && gConfig->lookupValue("rec.fit.y.p2", fFitPars[1][2])))
        Warning(here, "Cannot find setting \"rec.beamfit\", using default value ......");

    if (!(gConfig->lookupValue("rec.corr.t", fExtTgtCorrT)
            && gConfig->lookupValue("rec.corr.p", fExtTgtCorrP)
            && gConfig->lookupValue("rec.corr.d", fExtTgtCorrD)))
        Warning(here, "Cannot find setting \"rec.corr\", using default value ......");

    if (!gConfig->lookupValue("drift.llimit", fDriftLimit))
        Warning(here, "Cannot find setting \"drift.llimit\", using default value ......");

    if (fDebug > 0) {
        Info(here, "fHRSAngle\t= %le", fHRSAngle / kDEG);
        Info(here, "fSieveZ  \t= %le", fSieveZ);
    }

    return 0;
}

ClassImp(G2PRec)