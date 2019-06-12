#ifndef _MARS_MATERIAL_H
#define _MARS_MATERIAL_H

#include "material.h"

//////////////////////////////////////////////////////////
// CUSATIS/MARS MATERIAL

class MarsMaterial;
class MarsMaterialStatus : public DisMechMaterialStatus
{
private:
    double omega0, maxEpsT, maxEpsN, temp_maxEpsT, temp_maxEpsN;
    double damage, temp_damage;
    double Kt, Ks, L, nt;
    double RAND_H;

    double giveS0tension(double omega) const;
    double giveS0compression(double omega) const;
    void computeOmega0();
    void computeDamage(Vector strain);
    void computeKsAnsKt();
public:
    MarsMaterialStatus(MarsMaterial *m, Element *e);
    virtual ~MarsMaterialStatus() {};
    void init();
    virtual void update();
    virtual Vector giveSecantNormalShearStiffness() const;
    virtual Vector giveSecantNormalShearStiffness(const Vector &strain);
    virtual Vector giveStress(const Vector &strain);
    virtual double giveValue(string code) const;
};


class MarsMaterial : public DisMechMaterial
{
private:
    double ft, Gt;
    double fs, Gs, fc, Kc, beta, mu, nc;
    double Lcrs, Lcrt;
public:
    MarsMaterial() { name = "Mars material"; };
    ~MarsMaterial() {};
    void readFromLine(istringstream &iss);
    MaterialStatus *giveNewMaterialStatus(Element *e);
    double giveFt() { return ft; }
    double giveFs() { return fs; }
    double giveFc() { return fc; }
    double giveBeta() { return beta; }
    double giveMu() { return mu; }
    double giveNc() { return nc; }
    double giveLcrs() { return Lcrs; }
    double giveLcrt() { return Lcrt; }
    double giveKc() { return Kc; }

    virtual void init();
};


/*
 * //////////////////////////////////////////////////////////
 * //COUPLED DISCRETE MECHANICAL - TRANSPORT MATERIAL
 *
 * class DisMechLinMat;
 * class DisMechLinMatStatus: public MaterialStatus {
 * private:
 *
 * public:
 *  DisMechLinMatStatus();
 *  DisMechLinMatStatus(DisMechLinMat *m);
 *  ~DisMechLinMatStatus(){};
 * };
 *
 * class DisMechLinMat: public Material {
 * private:
 *  double E0, alpha, density;
 * public:
 *  DisMechLinMat();
 *  ~DisMechLinMat(){};
 *  void readFromLine(istringstream &iss);
 *  vector<double> giveStress(vector<double> strain) const;
 *  MaterialStatus* giveNewMaterialStatus();
 * };
 */

#endif /* _MARS_MATERIAL_H */
