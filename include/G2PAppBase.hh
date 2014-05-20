// -*- C++ -*-

/* class G2PAppBase
 * Abstract base class for g2p reconstruction tools.
 * It provides fundamental functions like parsing configuration files.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Add configure functions.
//   Feb 2014, C. Gu, Modified for G2PRec.
//

#ifndef G2P_APPBASE_H
#define G2P_APPBASE_H

#include "TObject.h"

class G2PAppBase : public TObject {
public:
    virtual ~G2PAppBase();

    int GetDebugLevel() const;
    void SetDebugLevel(int level);

protected:
    G2PAppBase(); // No instance allowed for this class

    virtual int Initialize() = 0;
    virtual void Clear(Option_t* /*option*/ = "") = 0;

    // Geometry utility functions
    virtual void TCS2HCS(double x_tr, double y_tr, double z_tr, double angle, double &x_lab, double &y_lab, double &z_lab);
    virtual void TCS2HCS(double t_tr, double p_tr, double angle, double &t_lab, double &p_lab);
    virtual void HCS2TCS(double x_lab, double y_lab, double z_lab, double angle, double &x_tr, double &y_tr, double &z_tr);
    virtual void HCS2TCS(double t_lab, double p_lab, double angle, double &t_tr, double &p_tr);
    virtual void Project(double x, double y, double z, double zout, double t, double p, double &xout, double &yout);

    virtual void TRCS2FCS(const double* V5_tr, double angle, double* V5_fp);
    virtual void FCS2TRCS(const double* V5_fp, double angle, double* V5_tr);
    virtual void TRCS2DCS(const double* V5_tr, double angle, double* V5_det);
    virtual void DCS2TRCS(const double* V5_det, double angle, double* V5_tr);
    virtual void FCS2DCS(const double* V5_fp, double angle, double* V5_det);
    virtual void DCS2FCS(const double* V5_det, double angle, double* V5_fp);

    virtual int Configure();

    // General status variables
    int fDebug;

private:
    ClassDef(G2PAppBase, 1)
};

inline int G2PAppBase::GetDebugLevel() const
{
    return fDebug;
}

inline void G2PAppBase::SetDebugLevel(int level)
{
    fDebug = level;
}

#endif
