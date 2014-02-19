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

    delete[] pDrift;
    pDrift = NULL;
}

int G2PRec::Process(const double* V5bpm_lab, const double* V5tp_tr, double* V5rec_tr)
{
    static const char* const here = "Process()";

    if (fDebug > 2) Info(here, " ");

    fV5tpmat_tr[0] = GetEffBPMX(V5bpm_lab, V5tp_tr);
    fV5tpmat_tr[1] = V5tp_tr[1];
    fV5tpmat_tr[2] = V5tp_tr[2];
    fV5tpmat_tr[3] = V5tp_tr[3];
    fV5tpmat_tr[4] = V5tp_tr[4];

#ifdef USE_BPMY
    double x, y, z;
    HCS2TCS(V5bpm_lab[0], V5bpm_lab[2], V5bpm_lab[4], fHRSAngle, x, y, z);
    fV5tpmat_tr[2] = y;
#endif

    if (fDebug > 1) {
        Info(here, "tpmat_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tpmat_tr[0], fV5tpmat_tr[1], fV5tpmat_tr[2], fV5tpmat_tr[3], fV5tpmat_tr[4]);
    }

    Project(fV5tpmat_tr[0], fV5tpmat_tr[2], 0.0, fSieveZ, fV5tpmat_tr[1], fV5tpmat_tr[3], fV5sieveproj_tr[0], fV5sieveproj_tr[2]);
    fV5sieveproj_tr[1] = fV5tpmat_tr[1];
    fV5sieveproj_tr[3] = fV5tpmat_tr[3];
    fV5sieveproj_tr[4] = fV5tpmat_tr[4];

    if (fDebug > 1) {
        Info(here, "sivproj_tr: %10.3e %10.3e %10.3e %10.3e %10.3e", fV5sieveproj_tr[0], fV5sieveproj_tr[1], fV5sieveproj_tr[2], fV5sieveproj_tr[3], fV5sieveproj_tr[4]);
    }

    pDrift->Drift(fV5sieveproj_tr, fSieveZ, fHRSMomentum, fHRSAngle, fRecZ, fV5tprec_tr);

    if (fDebug > 1) {
        Info(here, "tprec_tr  : %10.3e %10.3e %10.3e %10.3e %10.3e", fV5tprec_tr[0], fV5tprec_tr[1], fV5tprec_tr[2], fV5tprec_tr[3], fV5tprec_tr[4]);
    }

    for (int i = 0; i < 5; i++) V5rec_tr[i] = fV5tprec_tr[i];

    return 0;
}

int G2PRec::Initialize()
{
    //static const char* const here = "Initialize()";

    return 0;
}

void G2PRec::Clear()
{
    memset(fV5tpmat_tr, 0, sizeof (fV5tpmat_tr));
    memset(fV5sieveproj_tr, 0, sizeof (fV5sieveproj_tr));
    memset(fV5tprec_tr, 0, sizeof (fV5tprec_tr));
}

double G2PRec::GetEffBPMX(const double* V5bpm_lab, const double* V5tp_tr)
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
    HCS2TCS(V5bpm_lab[0], V5bpm_lab[2], V5bpm_lab[4], fHRSAngle, V3bpm_tr[0], V3bpm_tr[1], V3bpm_tr[2]);

    double xbpm_tr = V3bpm_tr[0];
    double xeffbpm_tr = xbpm_tr;

    if (fFieldRatio > 1e-5) {
        double p = (1 + V5tp_tr[4]) * fHRSMomentum;
        if (fFieldRatio < 0.75) xeffbpm_tr = xbpm_tr - (3.14345 / p + 0.0183611) / 1000 * fFieldRatio / 0.5;
        else xeffbpm_tr = xbpm_tr - (6.11766 / p + 0.14139) / 1000 * fFieldRatio / 1.0;
    }

    if (fDebug > 2) {
        Info(here, "%10.3e", xeffbpm_tr);
    }

    return xeffbpm_tr;
}

int G2PRec::Configure()
{
    static const char* const here = "Configure()";

    if (!gConfig->lookupValue("field.ratio", fFieldRatio))
        Warning(here, "Cannot find setting \"field.ratio\", using default value ......");

    if (!gConfig->lookupValue("hrs.angle", fHRSAngle))
        Warning(here, "Cannot find setting \"hrs.angle\", using default value ......");

    if (!gConfig->lookupValue("sieve.z", fSieveZ))
        Warning(here, "Cannot find setting \"sieve.z\", using default value ......");

    if (!gConfig->lookupValue("rec.z", fRecZ))
        Warning(here, "Cannot find setting \"rec.z\", using default value ......");

    return 0;
}

ClassImp(G2PRec)