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

#define USE_BPMY 1

using namespace std;
using namespace libconfig;

static const double kDEG = 3.14159265358979323846 / 180.0;

G2PRec::G2PRec() :
fHRSAngle(5.767 * kDEG), fHRSMomentum(2.254), fFieldRatio(0.0), fSieveZ(800.0), fRecZ(0.0)
{
    // Constructor

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

int G2PRec::Process(const float* V5bpm_bpm, const double* V5tp_tr, double * V5rec_tr, double* V5rec_lab)
{
    static const char* const here = "Process()";

    fV5bpm_bpm[0] = V5bpm_bpm[0] / 1000.0;
    fV5bpm_bpm[1] = V5bpm_bpm[1];
    fV5bpm_bpm[2] = V5bpm_bpm[2] / 1000.0;
    fV5bpm_bpm[3] = V5bpm_bpm[3];
    fV5bpm_bpm[4] = V5bpm_bpm[4] / 1000.0;

    fV5tpmat_tr[0] = V5tp_tr[0];
    fV5tpmat_tr[1] = V5tp_tr[1];
    fV5tpmat_tr[2] = V5tp_tr[2];
    fV5tpmat_tr[3] = V5tp_tr[3];
    fV5tpmat_tr[4] = V5tp_tr[4];

    //fV5tpmat_tr[0] = GetEffBPMX();

#ifdef USE_BPMY
    double V3temp_tr[3];
    HCS2TCS(fV5bpm_bpm[0], fV5bpm_bpm[2], fV5bpm_bpm[4], fHRSAngle, V3temp_tr[0], V3temp_tr[1], V3temp_tr[2]);
    fV5tpmat_tr[2] = V3temp_tr[1];
#endif

    if (fDebug > 0) {
        Info(here, "bpm_bpm   : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5bpm_bpm[0], fV5bpm_bpm[1], fV5bpm_bpm[2], fV5bpm_bpm[3], fV5bpm_bpm[4]);
        Info(here, "tpmat_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpmat_tr[0], fV5tpmat_tr[1], fV5tpmat_tr[2], fV5tpmat_tr[3], fV5tpmat_tr[4]);
    }

    Project(fV5tpmat_tr[0], fV5tpmat_tr[2], 0.0, fSieveZ, fV5tpmat_tr[1], fV5tpmat_tr[3], fV5sieveproj_tr[0], fV5sieveproj_tr[2]);
    fV5sieveproj_tr[1] = fV5tpmat_tr[1];
    fV5sieveproj_tr[3] = fV5tpmat_tr[3];
    fV5sieveproj_tr[4] = fV5tpmat_tr[4];

    if (fDebug > 0) {
        Info(here, "sivproj_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieveproj_tr[0], fV5sieveproj_tr[1], fV5sieveproj_tr[2], fV5sieveproj_tr[3], fV5sieveproj_tr[4]);
    }

    pDrift->Drift(fV5sieveproj_tr, fSieveZ, fHRSMomentum, fHRSAngle, 0.0, fV5tprec_tr);
    TCS2HCS(fV5tprec_tr[0], fV5tprec_tr[2], 0.0, fHRSAngle, fV5tprec_lab[0], fV5tprec_lab[2], fV5tprec_lab[4]);
    TCS2HCS(fV5tprec_tr[1], fV5tprec_tr[3], fHRSAngle, fV5tprec_lab[1], fV5tprec_lab[3]);

    double x[3] = {fV5tprec_lab[0], fV5tprec_lab[2], fV5tprec_lab[4]};
    double p[3] = {fHRSMomentum * (1 + fV5tprec_tr[4]) * sin(fV5tprec_lab[1]) * cos(fV5tprec_lab[3]),
        fHRSMomentum * (1 + fV5tprec_tr[4]) * sin(fV5tprec_lab[1]) * sin(fV5tprec_lab[3]),
        fHRSMomentum * (1 + fV5tprec_tr[4]) * cos(fV5tprec_lab[1])};
    pDrift->Drift(x, p, fRecZ, x, p);
    double temp;
    fV5tprec_lab[0] = x[0];
    fV5tprec_lab[1] = acos(p[2] / (fHRSMomentum * (1 + fV5tprec_tr[4])));
    fV5tprec_lab[2] = x[1];
    fV5tprec_lab[3] = atan2(p[1], p[0]);
    fV5tprec_lab[4] = x[2];
    HCS2TCS(x[0], x[1], x[2], fHRSAngle, fV5tprec_tr[0], fV5tprec_tr[2], temp);
    HCS2TCS(fV5tprec_lab[1], fV5tprec_lab[2], fHRSAngle, fV5tprec_tr[1], fV5tprec_tr[3]);

    if (fDebug > 0) {
        Info(here, "tprec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tprec_tr[0], fV5tprec_tr[1], fV5tprec_tr[2], fV5tprec_tr[3], fV5tprec_tr[4]);
        Info(here, "tprec_lab : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tprec_lab[0], fV5tprec_lab[1], fV5tprec_lab[2], fV5tprec_lab[3], fV5tprec_lab[4]);
    }

    for (int i = 0; i < 5; i++) {
        V5rec_tr[i] = fV5tprec_tr[i];
        V5rec_lab[i] = fV5tprec_lab[i];
    }

    return 0;
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

void G2PRec::Clear()
{
    memset(fV5bpm_bpm, 0, sizeof (fV5bpm_bpm));
    memset(fV5tpmat_tr, 0, sizeof (fV5tpmat_tr));
    memset(fV5sieveproj_tr, 0, sizeof (fV5sieveproj_tr));
    memset(fV5tprec_tr, 0, sizeof (fV5tprec_tr));
    memset(fV5tprec_lab, 0, sizeof (fV5tprec_lab));
}

double G2PRec::GetEffBPMX()
{
    static const char* const here = "GetEffBPMX()";

    // Fit result:
    //
    // (Xbpm_tr-Xtg_tr) vs Z
    // ([0]+[1]*x)
    // Fitting result of (Xbpm_tr-Xtg_tr) vs Z @ 2.5T
    // p0                        =    -0.011838   +/-   0.00132798
    // p1                        =      49.856    +/-   0.163115
    //
    // (Xtg_tr-Xtgproj_tr) vs P
    // ([0]+[1]/x)
    // Fitting result of (Xtg_tr-Xtgproj_tr) vs P @ 2.5T
    // p0                        =    0.0183611   +/-   0.0105237
    // p1                        =      3.14345   +/-   0.0105453
    // Fitting result of (Xtg_tr-Xtgproj_tr) vs P @ 5.0T
    // p0                        =      0.14139   +/-   0.018683
    // p1                        =      6.11766   +/-   0.0187211
    //

    double V3bpm_tr[5];
    HCS2TCS(fV5bpm_bpm[0], fV5bpm_bpm[2], fV5bpm_bpm[4], fHRSAngle, V3bpm_tr[0], V3bpm_tr[1], V3bpm_tr[2]);

    double xbpm_tr = V3bpm_tr[0];
    double xeffbpm_tr = xbpm_tr;

    if (fFieldRatio > 1e-5) {
        double p = (1 + fV5tpmat_tr[4]) * fHRSMomentum;
        if (fFieldRatio < 0.75) xeffbpm_tr = xbpm_tr - (3.14345 / p + 0.0183611) / 1000 * fFieldRatio / 0.5;
        else xeffbpm_tr = xbpm_tr - (6.11766 / p + 0.14139) / 1000 * fFieldRatio / 1.0;
    }

    if (fDebug > 1) {
        Info(here, "%10.3e", xeffbpm_tr);
    }

    return xeffbpm_tr;
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

    if (!gConfig->lookupValue("rec.z", fRecZ))
        Warning(here, "Cannot find setting \"rec.z\", using default value ......");

    if (fDebug > 0) {
        Info(here, "fHRSAngle\t= %le", fHRSAngle / kDEG);
        Info(here, "fSieveZ  \t= %le", fSieveZ);
        Info(here, "fRecZ    \t= %le", fRecZ);
    }

    return 0;
}

ClassImp(G2PRec)