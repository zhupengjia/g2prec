// -*- C++ -*-

/* class G2PField
 * Generate field map for HallB magnets.
 * Calculate field strength of a particular point from the field map using Lagrange polynomial interpolation, default is 2nd order.
 *
 * Field map may have an angle respect to the lab coordinates.
 * Use SetEulerAngle() to set this angle and the program will rotate the field map to correct direction.
 */

// History:
//   Mar 2013, C. Gu, First public version.
//   Sep 2013, C. Gu, Put HallB field map into G2PField.
//   Feb 2014, C. Gu, Modified for G2PRec.
//

#ifndef G2P_FIELD_H
#define G2P_FIELD_H

#include <vector>

#include "G2PAppBase.hh"

using namespace std;

class G2PField : public G2PAppBase {
public:
    G2PField();
    virtual ~G2PField();

    virtual void GetField(const double* x, double* b);

protected:
    virtual int Initialize();
    virtual void Clear();

    void SetRotationMatrix();

    virtual int ReadMap();
    virtual int CreateMap();

    virtual int Interpolate(const double* x, double* b, int order);

    void TransLab2Field(const double* x, double* xout);
    void TransField2Lab(const double* b, double* bout);

    virtual int Configure();

    const char* fMapFile;

    vector<vector<vector<double> > > fBField;

    double fOrigin[3];

    double fZMin, fZMax;
    double fRMin, fRMax;
    double fZStep, fRStep;
    int nZ, nR;

    bool fRotation;
    double fEulerAngle[3];
    double fRotationMatrix[2][3][3];

    double fRatio;

private:
    ClassDef(G2PField, 1)
};

#endif