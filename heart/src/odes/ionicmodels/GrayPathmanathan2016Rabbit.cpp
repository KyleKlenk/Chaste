/*

Copyright (c) 2005-2019, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/
#include "GrayPathmanathan2016Rabbit.hpp"
#include "OdeSystemInformation.hpp"
#include <cmath>

//
// Model-scope constant parameters
//
const double GrayPathmanathan2016Rabbit::mg_Na = 11.0 ;  // ms/cm^2
const double GrayPathmanathan2016Rabbit::mE_Na = 65.0;  // mV
const double GrayPathmanathan2016Rabbit::mE_K = -83.0;    //mV
const double GrayPathmanathan2016Rabbit::mE_h = -74.7 ;   // mV
const double GrayPathmanathan2016Rabbit::mE_m = -41.0; // (mV)
const double GrayPathmanathan2016Rabbit::mk_m = -4.0; //  (mV)  possible typo in Sebastian's code
const double GrayPathmanathan2016Rabbit::mk_r = 21.2766; // (mV) as in 1/b in GP2016 paper
const double GrayPathmanathan2016Rabbit::mk_h = 4.4; // (mV)
const double GrayPathmanathan2016Rabbit::mtau_m = 0.12; // (ms)
const double GrayPathmanathan2016Rabbit::mtau_ho = 6.80738; // (ms)
const double GrayPathmanathan2016Rabbit::mdelta_h = 0.799163; // (dimensionless)
const double GrayPathmanathan2016Rabbit::mg_K = 0.3; // (mS/cm^2)
const double GrayPathmanathan2016Rabbit::mC_m = 1.0; // (uF/cm^2)
//const double GrayPathmanathan2016Rabbit::mstim_start = 10.0; // (ms)
//const double GrayPathmanathan2016Rabbit::mstim_period = 1000.0; // (ms)
//const double GrayPathmanathan2016Rabbit::mstim_duration = 1.0; // (ms)
//const double GrayPathmanathan2016Rabbit::mstim_amplitude = 80.0; // (uA/uF)




GrayPathmanathan2016Rabbit::GrayPathmanathan2016Rabbit(
        boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
        boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus)
    : AbstractCardiacCell(pOdeSolver, 3, 0, pIntracellularStimulus)
{
    mpSystemInfo = OdeSystemInformation<GrayPathmanathan2016Rabbit>::Instance();

    Init();
}

GrayPathmanathan2016Rabbit::~GrayPathmanathan2016Rabbit(void)
{
}

void GrayPathmanathan2016Rabbit::EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY)
{
    double membrane_V = rY[0]; // v
    double variable_m = rY[1]; // m
    double variable_h = rY[2]; // h


    // Compute algebraic variables: 
    double m_inf = 1.0/(1.0+exp((membrane_V-mE_m)/(mk_m))); // dimensionless
    double h_inf = 1.0/(1.0+exp((membrane_V-mE_h)/(mk_h))); // dimensionless
    double tau_h = 2.0*mtau_ho*exp((mdelta_h*(membrane_V-mE_h))/mk_h)/(1+exp((membrane_V-mE_h)/(mk_h))); // ms
    double i_stim = GetIntracellularAreaStimulus(time);
  

    // dV/dt
    double membrane_V_prime = 0;
    // do not update voltage if the mSetVoltageDerivativeToZero flag has been set
    if (!mSetVoltageDerivativeToZero)
    {
    double i_na = mg_Na*pow(variable_m,3)*variable_h*(membrane_V - mE_Na); // uA/cm^2 
    double i_k = mg_K*(membrane_V-mE_K)*exp(-(membrane_V-mE_K)/mk_r); // uA/cm^2  
    membrane_V_prime = -(i_na+i_k)/mC_m+i_stim;

    }

    // dm/dt, dh/dt
    double variable_m_prime =  (m_inf-variable_m)/mtau_m;
    double variable_h_prime =  (h_inf-variable_h)/tau_h;
    
   // std::cout<< membrane_V_prime << std::endl;     

    rDY[0] = membrane_V_prime;
    rDY[1] = variable_m_prime;
    rDY[2] = variable_h_prime;
}

double GrayPathmanathan2016Rabbit::GetIIonic(const std::vector<double>* pStateVariables)
{
    if (!pStateVariables) pStateVariables = &mStateVariables;
    double membrane_V = (*pStateVariables)[mVoltageIndex];
    double variable_m = (*pStateVariables)[1];
    double variable_h = (*pStateVariables)[2];
    double i_na = mg_Na*pow(variable_m,3)*variable_h*(membrane_V - mE_Na); // uA/cm^2 
    double i_k = mg_K*(membrane_V-mE_K)*exp(-(membrane_V-mE_K)/mk_r); // uA/cm^2 
    double fake_ionic_current = (i_na+i_k)/mC_m;
    return fake_ionic_current;
}

template<>
void OdeSystemInformation<GrayPathmanathan2016Rabbit>::Initialise(void)
{
    /*
     * State variables
     */
    this->mVariableNames.push_back("V");
    this->mVariableUnits.push_back("mV");
    this->mInitialConditions.push_back(-83.0);

    this->mVariableNames.push_back("m");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(1.0/(1.0+exp((-83.0+41.0)/(-4.0))));

    this->mVariableNames.push_back("h");
    this->mVariableUnits.push_back("");
    this->mInitialConditions.push_back(1.0/(1.0+exp((-83.0+74.7)/(4.4))));

    this->mInitialised = true;
}