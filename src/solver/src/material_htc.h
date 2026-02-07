#ifndef _MATERIAL_HTC_H
#define _MATERIAL_HTC_H

#include "material_vectorial.h"

//////////////////////////////////////////////////////////
// HTC MATERIAL

class HTCMaterial;
class HTCMaterialStatus : public TensTrsprtMaterialStatus
{
private:
    double temp_alphac, temp_alphas, alphac, alphas;
    double Dh, qh, qT, dwe_dh;
    double h, T;

    void updateMaterialParameters(double timeStep);
public:
    HTCMaterialStatus(HTCMaterial *m, Element *e, unsigned ipnum);
    virtual ~HTCMaterialStatus() {};
    virtual void init();
    virtual void update();
    virtual Matrix giveStiffnessTensor(std :: string type) const;
    virtual void computeStress( double timeStep);
    virtual void computeStressWithFrozenIntVars( double timeStep);
    virtual bool giveValues(std :: string code, Vector &result) const;
    virtual Vector giveInternalSource() const;
    virtual Matrix giveDampingTensor() const;
    virtual Matrix giveMassTensor() const;
    virtual void setParameterValue(std :: string code, double value);
};


class HTCMaterial : public TensTrsprtMaterial
{
private:
    Matrix permeabilityTensor;
    double kappa, D1, rho, c, ct, Qcinf, Qsinf, s, EacR, EasR, Ac1, Ac2, As1, As2, alphacinf, alphasinf, a, b, etas, etac, kcvg, ksvg, g1, kappac, D0, EadR, T0, w0, n;
    double init_alphac, init_alphas;

public:
    HTCMaterial(unsigned dimension) : TensTrsprtMaterial(dimension) { name = "HTC material"; produceInternalSources = true; dampMatUpdate = true; };
    virtual ~HTCMaterial() {};
    virtual void readFromLine(std :: istringstream &iss);
    virtual MaterialStatus *giveNewMaterialStatus(Element *e, unsigned ipnum);
    virtual void init(MaterialContainer *matcont);
    Matrix givePermeabilityTensor() const { return permeabilityTensor; };
    double giveKappa()const { return kappa; };
    double giveD1()const { return D1; };
    double giveRho()const { return rho; };
    double giveC()const { return c; };
    double giveQcinf()const { return Qcinf; };
    double giveQsinf()const { return Qsinf; };
    double giveS()const { return s; };
    double giveEacR()const { return EacR; };
    double giveEasR()const { return EasR; };
    double giveAc1()const { return Ac1; };
    double giveAc2()const { return Ac2; };
    double giveAs1()const { return As1; };
    double giveAs2()const { return As2; };
    double giveAlphacinf()const { return alphacinf; };
    double giveAlphasinf()const { return alphasinf; };
    double giveA()const { return a; };
    double giveB()const { return b; };
    double giveEtas()const { return etas; };
    double giveEtac()const { return etac; };
    double giveKcvg()const { return kcvg; };
    double giveKsvg()const { return ksvg; };
    double giveG1()const { return g1; };
    double giveKappac()const { return kappac; };
    double giveD0()const { return D0; };
    double giveEadR()const { return EadR; };
    double giveT0()const { return T0; };
    double giveCt()const { return ct; };
    double giveN()const { return n; };
    double giveW0()const { return w0; };
    double giveInitAlphac()const { return init_alphac; };
    double giveInitAlphas()const { return init_alphas; };
};
#endif /* _MATERIAL_HTC_H */
