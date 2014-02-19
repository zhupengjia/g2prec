// -*- C++ -*-

/* class G2PDrift
 * Use Nystrom-Runge-Kutta method to derive the trajectory of a charged particle in static magnetic field.
 * Drift() has 2 prototypes, one for lab coordinates and one for HRS transportation coordinates.
 */

// History:
//   Feb 2013, C. Gu, First public version.
//   Feb 2013, C. Gu, Change algorithm to Nystrom-Runge-Kutta method.
//   Mar 2013, C. Gu, Add flexible step length and boundary check.
//   Oct 2013, J. Liu, Add drift function to stop at a cylinder boundary.
//   Feb 2014, C. Gu, Modified for G2PRec.
//

#ifndef G2P_DRIFT_H
#define G2P_DRIFT_H

#include "G2PAppBase.hh"

class G2PField;

class G2PDrift : public G2PAppBase {
public:
    G2PDrift();
    virtual ~G2PDrift();

    virtual double Drift(const double* x, const double* p, double zf, double *xout, double *pout); // HCS
    virtual double Drift(const double* x, double z_tr, double p, double angle, double zf_tr, double* xout); // TCS

protected:
    typedef double (G2PDrift::*pfDriftHCS_)(const double*, const double*, double, double*, double*);
    typedef double (G2PDrift::*pfDriftTCS_)(const double*, double, double, double, double, double*);

    virtual int Initialize();
    virtual void Clear();

    double DriftHCS(const double* x, const double* p, double zf, double *xout, double *pout);
    double DriftHCSNF(const double* x, const double* p, double zf, double *xout, double *pout);

    double DriftTCS(const double* x, double z_tr, double p, double angle, double zf_tr, double* xout);
    double DriftTCSNF(const double* x, double z_tr, double p, double angle, double zf_tr, double* xout);

    void NystromRK4(const double* x, const double* dxdt, double step, double* xo, double* err);
    double DistChord();
    void ComputeRHS(const double* x, double* dxdt);

    virtual int Configure();

    double fM0;
    double fQ, fQSave;
    double fStep, fStepLimit, fErrLoLimit, fErrHiLimit;
    double fVelocity, fVelocity2, fGamma;
    double fCof;
    double fField[3];
    double fIPoint[3];
    double fMPoint[3];
    double fEPoint[3];

    G2PField* pField;
    double fFieldRatio;

    pfDriftHCS_ pfDriftHCS;
    pfDriftTCS_ pfDriftTCS;

private:
    ClassDef(G2PDrift, 1)
};

#endif