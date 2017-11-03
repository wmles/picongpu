/* Copyright 2013-2017 Axel Huebl, Heiko Burau, Rene Widera, Richard Pausch, Stefan Tietze
 *
 * This file is part of PIConGPU.
 *
 * PIConGPU is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * PIConGPU is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with PIConGPU.
 * If not, see <http://www.gnu.org/licenses/>.
 */



#pragma once

#include <pmacc/types.hpp>
#include "picongpu/simulation_defines.hpp"

namespace picongpu
{
/** Wavepaket with spacial gaussian envelope and complicated longitudinal shape.
 * Allows prepulse, and also defining ASE with three (t, intensity)-points, where time is counted from the very beginning and the intensity is in multiples of the main peak.
 *
 *  _not_ including correction to laser formular derived from vector potential, so the integration
 *  along propagation direction gives 0
 *  this is (would be^^) important for few-cycle laser pulses
 */
namespace laserIlja
{
    // various computations 
    constexpr float_X laserTimeShift = laser::initPlaneY * CELL_HEIGHT / SPEED_OF_LIGHT;
    constexpr float_64 f = SPEED_OF_LIGHT / WAVE_LENGTH;
    constexpr float_64 w = 2.0 * PI * f;
    constexpr float_64 endUpramp = -0.5 * LASER_NOFOCUS_CONSTANT;
    constexpr float_64 startDownramp = 0.5 * LASER_NOFOCUS_CONSTANT;

    // helper functions, called from laserLongitudinal
    float_X gauss(float_64 t)
    {
        return math::exp(-0.25 * (t / PULSE_LENGTH) * (t / PULSE_LENGTH));
    }

    float_X extrapoliere(float_64 t1, float_X a1, float_64 t2, float_X a2, float_64 t)
    {
        float_X log1 = (t2 - t) * math::log(a1);
        float_X log2 = (t - t1) * math::log(a2);
        return math::exp((log1 + log2)/(t2 - t1));
    }

    float_X get_envelope(float_64 runTime)
    {
        float_X env = 0.0;
        if ((-0.5 * RAMP_INIT * PULSE_LENGTH < runTime) and (runTime < 0.))
            env = AMP_1 * gauss(runTime - 0.) + AMP_PREPULSE * gauss(runTime - TIME_PREPULSE);
        else if ((0. <= runTime) and (runTime < TIME_PEAKPULSE))
        {
            float_X ase_zum_peakpulse = extrapoliere(TIME_2, AMP_2, TIME_3, AMP_3, TIME_PEAKPULSE) / AMPLITUDE;
            if (ase_zum_peakpulse > 0.5) throw std::invalid_argument("\n\nAttention, the intensities of the ASE are very large, the extrapolation to the time of the main pulse would give more than 50% of the pulse amplitude - this is not an ASE physically, probably something wrong?!\n");
            env += AMPLITUDE * (1.-ase_zum_peakpulse) * gauss(runTime-TIME_PEAKPULSE);
            env += AMP_PREPULSE * gauss(runTime - TIME_PREPULSE);
            if ((TIME_1 < runTime) and (runTime < TIME_2))
                env += extrapoliere(TIME_1, AMP_1, TIME_2, AMP_2, runTime);
            else
                env += extrapoliere(TIME_2, AMP_2, TIME_3, AMP_3, runTime);
        }
        else if (TIME_PEAKPULSE <= runTime) env = AMPLITUDE * gauss(runTime-TIME_PEAKPULSE);
        return env;
    }

    HINLINE float3_X laserLongitudinal(uint32_t currentStep, float_X& phase)
    {
        float_X envelope;
        float3_X elong(float3_X::create(0.0));

        // a symmetric pulse will be initialized at position z=0 for
        // a time of RAMP_INIT * PULSE_LENGTH + LASER_NOFOCUS_CONSTANT = INIT_TIME.
        // we shift the complete pulse for the half of this time to start with
        // the front of the laser pulse.

        /* initialize the laser not in the first cell is equal to a negative shift
         * in time
         */
        const float_64 runTime = (DELTA_T * currentStep - laserTimeShift - 0.5 * RAMP_INIT * PULSE_LENGTH);

        const float_64 tau = PULSE_LENGTH * sqrt(2.0);

        phase += float_X(w * runTime) + LASER_PHASE ;

        envelope = get_envelope(runTime);


        float_64 correctionFactor = 0.0;
        /* I don't want to delete the correctionFactor and the old implementation using it.
         * Would be nice to use it, at least in the regions where gauss is dominating.
         */

        /* if (runTime > startDownramp)
        {
            // downramp = end
            const float_64 exponent =
                ((runTime - startDownramp)
                 / PULSE_LENGTH / sqrt(2.0));
            envelope *= math::exp(-0.5 * exponent * exponent);
            correctionFactor = (runTime - startDownramp)/(tau*tau*w);
        }
        else if(runTime < endUpramp)
        {
            // upramp = start
            const float_64 exponent = ((runTime - endUpramp) / PULSE_LENGTH / sqrt(2.0));
            envelope *= math::exp(-0.5 * exponent * exponent);
            correctionFactor = (runTime - endUpramp)/(tau*tau*w);
        } */

        if( Polarisation == LINEAR_X )
        {
            elong.x() = float_X(envelope * (math::sin(phase) + correctionFactor * math::cos(phase)));
        }
        else if( Polarisation == LINEAR_Z )
        {
            elong.z() = float_X(envelope * (math::sin(phase) + correctionFactor * math::cos(phase)));
        }
        else if( Polarisation == CIRCULAR )
        {
            elong.x() = float_X(envelope / sqrt(2.0) * (math::sin(phase) + correctionFactor * math::cos(phase)));
            elong.z() = float_X(envelope / sqrt(2.0) * (math::cos(phase) + correctionFactor * math::sin(phase)));
        }

        return elong;
    }

    /**
     *
     * @param elong
     * @param phase
     * @param posX
     * @param posZ
     * @return
     */
    HDINLINE float3_X laserTransversal(float3_X elong, const float_X, const float_X posX, const float_X posZ)
    {

        const float_X exp_x = posX * posX / (W0_X * W0_X);
        const float_X exp_z = posZ * posZ / (W0_Z * W0_Z);

        return elong * math::exp( float_X( -1.0 ) * ( exp_x + exp_z ) );

    }

}
}




