/*
 * Copyright (c) 2020-2024 Key4hep-Project.
 *
 * This file is part of Key4hep.
 * See https://key4hep.github.io/key4hep-doc/ for further info.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */


#include "GenfitMaterialInterface.h"

namespace GenfitInterface{

    double MeanExcEnergy_get(int Z);
    double MeanExcEnergy_get(TGeoMaterial*);
    GenfitMaterialInterface* GenfitMaterialInterface::m_instance = nullptr;

    GenfitMaterialInterface::GenfitMaterialInterface(
            const dd4hep::Detector* dd4hepGeo)
    {
        assert(nullptr!=dd4hepGeo);
        m_geoManager=&(dd4hepGeo->manager());
        assert(nullptr!=m_geoManager);
        debugLvl_=0;
    }

    GenfitMaterialInterface* GenfitMaterialInterface::getInstance(
            const dd4hep::Detector* dd4hepGeo)
    {
        if(nullptr == m_instance){
            m_instance = new GenfitMaterialInterface(dd4hepGeo);
        }
        genfit::MaterialEffects::getInstance()->init(m_instance);
        return m_instance;
    }

    void GenfitMaterialInterface::destruct(){
        delete m_instance;
    }

    TGeoManager* GenfitMaterialInterface::getGeoManager()
    {
        return m_geoManager;
    }

    bool GenfitMaterialInterface::initTrack(double posX, double posY,double posZ, 
                                            double dirX, double dirY, double dirZ)
    {
        // Move to the new point.
        bool result = !isSameLocation(posX, posY, posZ, kTRUE);

        // Set the intended direction.
        setCurrentDirection(dirX, dirY, dirZ);

        if (debugLvl_ > 0) {

            std::cout << "      GenfitMaterialInterface::initTrack at \n";
            std::cout << "      position:  "; TVector3(posX, posY, posZ).Print();
            std::cout << "      direction: "; TVector3(dirX, dirY, dirZ).Print();

        }

        return result;
    }


    genfit::Material GenfitMaterialInterface::getMaterialParameters()
    {
        TGeoMaterial* mat = getGeoManager()->GetCurrentVolume()->GetMedium()->GetMaterial();
    
        return genfit::Material(mat->GetDensity(), mat->GetZ(), mat->GetA(),
               mat->GetRadLen(), MeanExcEnergy_get(mat));
        
    }

    double GenfitMaterialInterface::findNextBoundary(   const genfit::RKTrackRep* rep,
                                                        const genfit::M1x7& stateOrig,
                                                        double sMax, // signed
                                                        bool varField)
    {
        
        const double delta(1.E-2); // cm, distance limit beneath which straight-line steps are taken.
        const double epsilon(1.E-1); // cm, allowed upper bound on arch
        // deviation from straight line

        genfit::M1x3 SA;
        genfit::M1x7 state7, oldState7;
        oldState7 = stateOrig;

        int stepSign(sMax < 0 ? -1 : 1);

        double s = 0;  // trajectory length to boundary

        const unsigned maxIt = 300;
        unsigned it = 0;

        // Initialize the geometry to the current location (set by caller).
        findNextBoundary(fabs(sMax) - s);
        double safety = getSafeDistance(); // >= 0
        double slDist = getStep();
    
        // this should not happen, but it happens sometimes.
        // The reason are probably overlaps in the geometry.
        // Without this "hack" many small steps with size of minStep will be made,
        // until eventually the maximum number of iterations are exceeded
        // and the extrapolation fails
        if (slDist < safety) slDist = safety;
        double step = slDist;

        if (debugLvl_ > 0)
            std::cout << "   safety = " << safety << "; step, slDist = " << step << "\n";

        while (1) {
            if (++it > maxIt){
                genfit::Exception exc("GenfitMaterialInterface::findNextBoundary ==> maximum number of iterations exceeded",__LINE__,__FILE__);
                exc.setFatal();
                throw exc;
            }

            // No boundary in sight?
            if (s + safety > fabs(sMax)) {
                if (debugLvl_ > 0)
                    std::cout << "   next boundary is further away than sMax \n";
                return stepSign*(s + safety); //sMax;
            }

            // Are we at the boundary?
            if (slDist < delta) {
                if (debugLvl_ > 0)
                    std::cout << "   very close to the boundary -> return"
                        << " stepSign*(s + slDist) = "
                        << stepSign << "*(" << s + slDist << ")\n";
                return stepSign*(s + slDist);
            }

            // We have to find whether there's any boundary on our path.

            // Follow curved arch, then see if we may have missed a boundary.
            // Always propagate complete way from original start to avoid
            // inconsistent extrapolations.
            state7 = stateOrig;
            rep->RKPropagate(state7, nullptr, SA, stepSign*(s + step), varField);

            if(debugLvl_){
                std::cout<<" RKTrackRep at state "<<state7[0]<<" "<<state7[1]
                    <<" "<<state7[2] <<" "<<state7[3] <<" "<<state7[4]
                    <<" "<<state7[5] <<" "<<state7[6]<<std::endl;
            }
            // Straight line distance² between extrapolation finish and
            // the end of the previously determined safe segment.
            double dist2 = (pow(state7[0] - oldState7[0], 2)
                    + pow(state7[1] - oldState7[1], 2)
                    + pow(state7[2] - oldState7[2], 2));
            // Maximal lateral deviation².
            double maxDeviation2 = 0.25*(step*step - dist2);

            if (step > safety
                    && maxDeviation2 > epsilon*epsilon) {
                // Need to take a shorter step to reliably estimate material,
                // but only if we didn't move by safety.  In order to avoid
                // issues with roundoff we check 'step > safety' instead of
                // 'step != safety'.  If we moved by safety, there couldn't have
                // been a boundary that we missed along the path, but we could
                // be on the boundary.

                // Take a shorter step, but never shorter than safety.
                step = std::max(step / 2, safety);
            } else {
                getGeoManager()->PushPoint();
                bool volChanged = initTrack(state7[0], state7[1], state7[2],
                        stepSign*state7[3], stepSign*state7[4],
                        stepSign*state7[5]);

                if (volChanged) {
                    if (debugLvl_ > 0)
                        std::cout << "   volChanged\n";
                    // Move back to start.
                    getGeoManager()->PopPoint();

                    // Extrapolation may not take the exact step length we asked
                    // for, so it can happen that a requested step < safety takes
                    // us across the boundary.  This is then the best estimate we
                    // can get of the distance to the boundary with the stepper.
                    if (step <= safety) {
                        if (debugLvl_ > 0)
                            std::cout <<
                                "   step <= safety, return stepSign*(s + step) = "
                                << stepSign*(s + step) << "\n";
                        return stepSign*(s + step);
                    }

                    // Volume changed during the extrapolation.  Take a shorter
                    // step, but never shorter than safety.
                    step = std::max(step / 2, safety);
                } else {
                    // we're in the new place, the step was safe, advance
                    s += step;

                    oldState7 = state7;
                    getGeoManager()->PopDummy();  // Pop stack, but stay in place.

                    findNextBoundary(fabs(sMax) - s);
                    safety = getSafeDistance();
                    slDist = getStep();

                    // this should not happen, but it happens sometimes.
                    // The reason are probably overlaps in the geometry.
                    // Without this "hack" many small steps with size of minStep will be made,
                    // until eventually the maximum number of iterations are exceeded and the extrapolation fails
                    if (slDist < safety)
                        slDist = safety;
                    step = slDist;
                    if (debugLvl_ > 0){
                        std::cout << "   safety = " << safety 
                            << "; step, slDist = " << step << "\n";
                    }
                }
            }
        }
    }


    /*
    Reference for elemental mean excitation energies at:
    http://physics.nist.gov/PhysRefData/XrayMassCoef/tab1.html

    Code ported from GEANT 3
    */

    const int MeanExcEnergy_NELEMENTS = 93; // 0 = vacuum, 1 = hydrogen, 92 = uranium
    const double MeanExcEnergy_vals[MeanExcEnergy_NELEMENTS] = {
        1.e15, 19.2, 41.8, 40.0, 63.7, 76.0, 78., 82.0,
        95.0, 115.0, 137.0, 149.0, 156.0, 166.0, 173.0, 173.0,
        180.0, 174.0, 188.0, 190.0, 191.0, 216.0, 233.0, 245.0,
        257.0, 272.0, 286.0, 297.0, 311.0, 322.0, 330.0, 334.0,
        350.0, 347.0, 348.0, 343.0, 352.0, 363.0, 366.0, 379.0,
        393.0, 417.0, 424.0, 428.0, 441.0, 449.0, 470.0, 470.0,
        469.0, 488.0, 488.0, 487.0, 485.0, 491.0, 482.0, 488.0,
        491.0, 501.0, 523.0, 535.0, 546.0, 560.0, 574.0, 580.0,
        591.0, 614.0, 628.0, 650.0, 658.0, 674.0, 684.0, 694.0,
        705.0, 718.0, 727.0, 736.0, 746.0, 757.0, 790.0, 790.0,
        800.0, 810.0, 823.0, 823.0, 830.0, 825.0, 794.0, 827.0,
        826.0, 841.0, 847.0, 878.0, 890.0, };
    // Logarithms of the previous table, used to calculate mixtures.
    const double MeanExcEnergy_logs[MeanExcEnergy_NELEMENTS] = {
        34.5388, 2.95491, 3.7329, 3.68888, 4.15418, 4.33073, 4.35671, 4.40672, 
        4.55388, 4.74493, 4.91998, 5.00395, 5.04986, 5.11199, 5.15329, 5.15329, 
        5.19296, 5.15906, 5.23644, 5.24702, 5.25227, 5.37528, 5.45104, 5.50126, 
        5.54908, 5.6058, 5.65599, 5.69373, 5.73979, 5.77455, 5.79909, 5.81114, 
        5.85793, 5.84932, 5.8522, 5.83773, 5.86363, 5.8944, 5.90263, 5.93754, 
        5.97381, 6.03309, 6.04973, 6.05912, 6.08904, 6.10702, 6.15273, 6.15273, 
        6.1506, 6.19032, 6.19032, 6.18826, 6.18415, 6.19644, 6.17794, 6.19032, 
        6.19644, 6.21661, 6.25958, 6.28227, 6.30262, 6.32794, 6.35263, 6.36303, 
        6.38182, 6.41999, 6.44254, 6.47697, 6.4892, 6.51323, 6.52796, 6.54247, 
        6.5582, 6.57647, 6.58893, 6.60123, 6.61473, 6.62936, 6.67203, 6.67203, 
        6.68461, 6.69703, 6.71296, 6.71296, 6.72143, 6.71538, 6.67708, 6.7178, 
        6.71659, 6.73459, 6.7417, 6.77765, 6.79122, };


    double MeanExcEnergy_get(int Z, bool logs = false) {
        assert(Z >= 0 && Z < MeanExcEnergy_NELEMENTS);
        if (logs)
            return MeanExcEnergy_logs[Z];
        else
            return MeanExcEnergy_vals[Z];
    }


    double MeanExcEnergy_get(TGeoMaterial* mat) {
        if (mat->IsMixture()) {
            // The mean excitation energy is calculated as the geometric mean
            // of the mEEs of the components, weighted for abundance.
            double logMEE = 0.;
            double denom  = 0.;
            TGeoMixture* mix = (TGeoMixture*)mat;
            for (int i = 0; i < mix->GetNelements(); ++i) {
                int index = int(mix->GetZmixt()[i]);
                double weight = mix->GetWmixt()[i] * mix->GetZmixt()[i] / mix->GetAmixt()[i];
                logMEE += weight * MeanExcEnergy_get(index, true);
                denom  += weight;
            }
            logMEE /= denom;
            return exp(logMEE);
        }

        // not a mixture
        int index = int(mat->GetZ());
        return MeanExcEnergy_get(index, false);
    }

    double GenfitMaterialInterface::getSafeDistance()
    {
        return getGeoManager()->GetSafeDistance();
    }

    double GenfitMaterialInterface::getStep()
    {
        return getGeoManager()->GetStep();
    }

    TGeoNode* GenfitMaterialInterface::findNextBoundary(double stepmax,
            const char* path,bool frombdr)
    {
        return getGeoManager()->FindNextBoundary(stepmax,path,frombdr);
    }

    bool GenfitMaterialInterface::isSameLocation(double posX, double posY,
            double posZ, bool change)
    {

        return getGeoManager()->IsSameLocation(
            posX,
            posY,
            posZ,change);
    }

    void GenfitMaterialInterface::setCurrentDirection(double nx, double ny,
            double nz)
    {
        return getGeoManager()->SetCurrentDirection(nx, ny, nz);

    }

}