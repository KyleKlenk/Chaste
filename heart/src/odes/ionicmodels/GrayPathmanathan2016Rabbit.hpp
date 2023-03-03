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


Siqi Wei
June 2022
Create cell model according to GrayPathmanathan2016Rabbit

There are a total of 17 entries in the constant variable array.

* VOI is time in component membrane (millisecond)
* STATES[0] is V in component membrane (mV)
* CONSTANTS[0] is g_Na (mS/mm^2)
* CONSTANTS[1] is E_Na (mV)
* CONSTANTS[2] is E_K (mV)
* CONSTANTS[3] is E_h (mV)
* CONSTANTS[4] is E_m (mV)
* CONSTANTS[5] is k_m (mV)
* CONSTANTS[6] is k_r (mV)
* CONSTANTS[7] is k_h (mV)
* CONSTANTS[8] is tau_m (ms)
* CONSTANTS[9] is tau_ho (ms)
* CONSTANTS[10] is delta_h (dimensionless)
* CONSTANTS[11] is g_K (mS/mm^2)
* CONSTANTS[12] is C_m (uF/mm^2)
* CONSTANTS[13] is stim_start (ms)
* CONSTANTS[14] is stim_period (ms)
* CONSTANTS[15] is stim_duration (ms)
* CONSTANTS[16] is stim_amplitude (uA/uF)
* STATES[1] is m (dimensionless)
* STATES[2] is h (dimensionless)
* RATES[0] is d/dt V (mV/ms)
* RATES[1] is d/dt m (1/ms)
* RATES[2] is d/dt h (1/ms)
* ALGEBRAIC[0] is m_inf (dimensionless)
* ALGEBRAIC[1] is h_inf (dimensionless)
* ALGEBRAIC[2] is tau_h (ms)
* ALGEBRAIC[3] is i_na (uA/mm^2)
* ALGEBRAIC[4] is i_k (uA/mm^2)
* ALGEBRAIC[5] is i_tot (uA/mm^2)
* ALGEBRAIC[6] is i_Stim (uA/uF)
*/

#ifndef _GRAYPATHMANATHAN2016RABBIT_HPP_
#define _GRAYPATHMANATHAN2016RABBIT_HPP_

#include "AbstractCardiacCell.hpp"
#include "AbstractStimulusFunction.hpp"
#include <vector>

/**
 * Represents the Gray-Pathmanathan2016 system of ODEs.
 */
class GrayPathmanathan2016Rabbit : public AbstractCardiacCell
{
private:
    /** constant parameters **/ 
    static const double mg_Na; // ms/cm^2 
    static const double mE_Na; // mV 
    static const double mE_K; //mV
    static const double mE_h; // mV
    static const double mE_m; // (mV)
    static const double mk_m; //  (mV)
    static const double mk_r; // (mV)
    static const double mk_h; // (mV)
    static const double mtau_m; // (ms)
    static const double mtau_ho; // (ms)
    static const double mdelta_h; // (dimensionless)
    static const double mg_K; // (mS/cm^2)
    static const double mC_m; // (uF/cm^2)
 //   static const double mstim_start; // (ms)
 //   static const double mstim_period; // (ms)
 //   static const double mstim_duration; // (ms)
 //  static const double mstim_amplitude; // (uA/uF)


public:
    /**
     * Constructor
     *
     * @param pOdeSolver is a pointer to the ODE solver
     * @param pIntracellularStimulus is a pointer to the intracellular stimulus
     */
    GrayPathmanathan2016Rabbit(boost::shared_ptr<AbstractIvpOdeSolver> pOdeSolver,
                                boost::shared_ptr<AbstractStimulusFunction> pIntracellularStimulus);

    /**
     * Destructor
     */
    ~GrayPathmanathan2016Rabbit();

    /**
     * Compute the RHS of the Gray-Pathmanathan2016 system of ODEs
     *
     * @param time  the current time, in milliseconds
     * @param rY  current values of the state variables
     * @param rDY  to be filled in with derivatives
     */
    void EvaluateYDerivatives(double time, const std::vector<double> &rY, std::vector<double>& rDY);

    /**
     * Calculates the ionic current.
     *
     * @param pStateVariables  optionally can be supplied to evaluate the ionic current at the
     *     given state; by default the cell's internal state will be used.
     *
     * @return the total ionic current
     */
    double GetIIonic(const std::vector<double>* pStateVariables=NULL);
};

#endif //_GRAYPATHMANATHAN2016RABBIT_HPP_
