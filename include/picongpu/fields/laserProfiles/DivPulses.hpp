/* Copyright 2018-2021 Ilja Goethel, Axel Huebl
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

#include "picongpu/simulation_defines.hpp"

#include <pmacc/mappings/simulation/SubGrid.hpp>
#include <pmacc/dataManagement/DataConnector.hpp>

#include "pmacc/misc/splitString.hpp"

namespace picongpu
{
    namespace fields
    {
        namespace laserProfiles
        {
            namespace divPulses
            {
                template<typename T_Params>
                struct Unitless : public T_Params
                {
                    using Params = T_Params;

                    static constexpr float_X WAVE_LENGTH
                        = float_X(Params::WAVE_LENGTH_SI / UNIT_LENGTH); // unit: meter
                    static constexpr float_X PULSE_LENGTH
                        = float_X(Params::PULSE_LENGTH_SI / UNIT_TIME); // unit: seconds (1 sigma)
                    static constexpr float_X LASER_NOFOCUS_CONSTANT
                        = float_X(Params::LASER_NOFOCUS_CONSTANT_SI / UNIT_TIME); // unit: seconds
                    static constexpr float_X AMPLITUDE
                        = float_X(Params::AMPLITUDE_SI / UNIT_EFIELD); // unit: Volt /meter
                    static constexpr float_X W0_X = float_X(Params::W0_X_SI / UNIT_LENGTH); // unit: meter
                    static constexpr float_X W0_Z = float_X(Params::W0_Z_SI / UNIT_LENGTH); // unit: meter

                    static constexpr float_64 TIME_START = float_64(Params::TIME_START_fs * 1.e-15 / UNIT_TIME); // start can be before TIME_0; start cuts away not only ramp but also pulses
                    static constexpr float_64 TIME_0 = float_64(Params::TIME_POINT_0_fs * 1.e-15 / UNIT_TIME);
                    static constexpr float_64 TIME_1 = float_64(Params::TIME_POINT_1_fs * 1.e-15 / UNIT_TIME);
                    static constexpr float_64 TIME_2 = float_64(Params::TIME_POINT_2_fs * 1.e-15 / UNIT_TIME);
                    static constexpr float_64 TIME_3 = float_64(Params::TIME_POINT_3_fs * 1.e-15 / UNIT_TIME);
                    static constexpr float_64 TIME_END = float_64(Params::TIME_ENDLASER_fs * 1.e-15 / UNIT_TIME);
                    static constexpr float_X startDownramp = LASER_NOFOCUS_CONSTANT;

                    // compile-time checks for physical sanity:
                    static_assert(
                        (TIME_0 <= TIME_1) && (TIME_1 <= TIME_2) && (TIME_2 <= TIME_3) && (TIME_3 <= 0.),
                        "The times in the parameters TIME_POINT_0/1/2/3 and the beginning of the plateau (which is at "
                        "time 0) should be in ascending order");

                   /* initialize the laser not in the first cell is equal to a negative shift
                     * in time
                     */
                    static constexpr float_X laserTimeShift = (Params::initPlaneY * CELL_HEIGHT / SPEED_OF_LIGHT);
                    static constexpr float_X INIT_TIME = float_X(TIME_END - TIME_START);

                    /* a symmetric pulse will be initialized at position z=0 for
                     * a time of RAMP_INIT * PULSE_LENGTH + LASER_NOFOCUS_CONSTANT = INIT_TIME.
                     * we shift the complete pulse for the half of this time to start with
                     * the front of the laser pulse.
                     */
                    static constexpr float_X time_start_init = TIME_START;
                    static constexpr float_64 f = SPEED_OF_LIGHT / WAVE_LENGTH;
                    static constexpr float_64 w = 2.0 * PI * f;
                };
            } // namespace divPulses

            namespace acc
            {
                template<typename T_Unitless>
                struct DivPulses : public T_Unitless
                {
                    using Unitless = T_Unitless;

                    float3_X m_elong;
                    float_X m_phase;
                    typename FieldE::DataBoxType m_dataBoxE;
                    DataSpace<simDim> m_offsetToTotalDomain;
                    DataSpace<simDim> m_superCellToLocalOriginCellOffset;

                    /** Device-Side Constructor
                     *
                     * @param superCellToLocalOriginCellOffset local offset in cells to current supercell
                     * @param offsetToTotalDomain offset to origin of global (@todo: total) coordinate system (possibly
                     * after transform to centered origin)
                     */
                    HDINLINE DivPulses(
                        typename FieldE::DataBoxType const& dataBoxE,
                        DataSpace<simDim> const& superCellToLocalOriginCellOffset,
                        DataSpace<simDim> const& offsetToTotalDomain,
                        float3_X const& elong)
                        : m_elong(elong)
                        , m_dataBoxE(dataBoxE)
                        , m_offsetToTotalDomain(offsetToTotalDomain)
                        , m_superCellToLocalOriginCellOffset(superCellToLocalOriginCellOffset)
                    {
                    }

                    /** device side manipulation for init plane (transversal)
                     *
                     * @tparam T_Args type of the arguments passed to the user manipulator functor
                     *
                     * @param cellIndexInSuperCell ND cell index in current supercell
                     */
                    template<typename T_Acc>
                    HDINLINE void operator()(T_Acc const&, DataSpace<simDim> const& cellIndexInSuperCell)
                    {
                        // coordinate system to global simulation as origin
                        DataSpace<simDim> const localCell(cellIndexInSuperCell + m_superCellToLocalOriginCellOffset);

                        // transform coordinate system to center of x-z plane of initialization
                        constexpr uint8_t planeNormalDir = 1u;
                        DataSpace<simDim> offsetToCenterOfPlane(m_offsetToTotalDomain);
                        offsetToCenterOfPlane[planeNormalDir] = 0; // do not shift origin of plane normal
                        floatD_X const pos
                            = precisionCast<float_X>(localCell + offsetToCenterOfPlane) * cellSize.shrink<simDim>();
                        // @todo add half-cells via traits::FieldPosition< Solver::NumicalCellType, FieldE >()

                        // transversal position only
                        float3_X const w0_3D(Unitless::W0_X, 0., Unitless::W0_Z);
                        auto const w0(w0_3D.shrink<simDim>().remove<planeNormalDir>());
                        auto const pos_trans(pos.remove<planeNormalDir>());
                        auto const exp_compos(pos_trans * pos_trans / (w0 * w0));
                        float_X const exp_arg(exp_compos.sumOfComponents());

                        m_elong *= math::exp(-1.0_X * exp_arg);

                        if(Unitless::initPlaneY != 0) // compile time if
                        {
                            /* If the laser is not initialized in the first cell we emit a
                             * negatively and positively propagating wave. Therefore we need to multiply the
                             * amplitude with a correction factor depending of the cell size in
                             * propagation direction.
                             * The negatively propagating wave is damped by the absorber.
                             *
                             * The `correctionFactor` assume that the wave is moving in y direction.
                             */
                            auto const correctionFactor = (SPEED_OF_LIGHT * DELTA_T) / CELL_HEIGHT * 2._X;

                            // jump over the guard of the electric field
                            m_dataBoxE(localCell + SuperCellSize::toRT() * GuardSize::toRT())
                                += correctionFactor * m_elong;
                        }
                        else
                        {
                            // jump over the guard of the electric field
                            m_dataBoxE(localCell + SuperCellSize::toRT() * GuardSize::toRT()) = m_elong;
                        }
                    }
                };
            } // namespace acc

            template<typename T_Params>
            struct DivPulses : public divPulses::Unitless<T_Params>
            {
                using Unitless = divPulses::Unitless<T_Params>;

                float3_X elong;
                float_X phase;
                typename FieldE::DataBoxType dataBoxE;
                DataSpace<simDim> offsetToTotalDomain;


                HDINLINE float_X gauss(
                    float_X const t,
                    float_X const rel_int,
                    float_X const t_center,
                    float_X const rel_len
                )
                {
                    float_X const exponent = (t - t_center) / float_X( rel_len * Unitless::PULSE_LENGTH );
                    return rel_int * math::exp(-0.5_X * exponent * exponent);
                }

                HDINLINE float_X sinpulse(
                    float_X const t,
                    float_X const rel_int,
                    float_X const t_start,
                    float_X const t_end
                )
                {
                    float_X const phase = M_PI * (t - t_start) / (t_end - t_start);
                    if ( phase < 0. || M_PI < phase )
                        return 0.;
                    else
                        return rel_int * math::sin(phase) * math::sin(phase);
                }

                /** get value of exponential curve through two points at given t
                 * t/t1/t2 given as float_X, since the envelope doesn't need the accuracy
                 */
                HDINLINE float_X extrapolate_expo(
                    float_X const t1,
                    float_X const a1,
                    float_X const t2,
                    float_X const a2,
                    float_X const t)
                {
                    const float_X log1 = (t2 - t) * math::log(a1);
                    const float_X log2 = (t - t1) * math::log(a2);
                    return math::exp((log1 + log2) / (t2 - t1));
                }

                HINLINE float_X get_envelope(float_X runTime)
                {
                    /* workaround for clang 5 linker issues
                     * `undefined reference to
                     * `picongpu::fields::laserProfiles::DivPulsesParam::INT_RATIO_POINT_1'`
                     */
                    float_X const INT_0 = Unitless::INT_RATIO_POINT_0;
                    float_X const INT_1 = Unitless::INT_RATIO_POINT_1;
                    float_X const INT_2 = Unitless::INT_RATIO_POINT_2;
                    float_X const INT_3 = Unitless::INT_RATIO_POINT_3;

                    float_X intens = 0.0;
                    bool const before_init = runTime < Unitless::time_start_init;
                    bool const before_ramp1 = runTime < Unitless::TIME_1;
                    bool const before_ramp2 = runTime < Unitless::TIME_2;
                    bool const before_plateau = runTime < 0.;
                    bool const after_plateau = Unitless::startDownramp <= runTime;
                    bool const after_endlaser = Unitless::TIME_END <= runTime;
                    bool const during_preramp = ( ! before_init ) && before_ramp1;
                    bool const during_ramp1 = ( ! before_ramp1 ) && before_ramp2;
                    bool const during_ramp2 = ( ! before_ramp2 ) && before_plateau;
                    bool const during_plateau = ( ! after_plateau ) && ( ! before_plateau );
                    bool const during_downramp = ( after_plateau ) && ( ! after_endlaser );
                    float_X const ramp_when_peakpulse
                        = extrapolate_expo(Unitless::TIME_2, INT_2, Unitless::TIME_3, INT_3, 0.);
                    if(ramp_when_peakpulse > 0.5)
                    {
                        log<picLog::PHYSICS>(
                            "Attention, the intensities of the laser upramp are very large! "
                            "The extrapolation of the last exponential to the time of "
                            "the peakpulse gives more than half of the amplitude of "
                            "the peak Gaussian. This is not a Gaussian at all anymore, "
                            "and physically very unplausible, check the params for misunderstandings!");
                    }

                    if(before_init)
                        intens = 0.;
                    else if(during_preramp)
                    {
                        intens = extrapolate_expo(Unitless::TIME_0, INT_0, Unitless::TIME_1, INT_1, runTime);
                    }
                    else if(during_ramp1)
                    {
                        intens = extrapolate_expo(Unitless::TIME_1, INT_1, Unitless::TIME_2, INT_2, runTime);
                    }
                    else if(during_ramp2)
                    {
                        intens = extrapolate_expo(Unitless::TIME_2, INT_2, Unitless::TIME_3, INT_3, runTime);
                    }

                    if(before_plateau)
                    {
                        intens += gauss(runTime, 1._X - ramp_when_peakpulse, 0., 1.);
                    }
                    else if(during_plateau)
                    {
                        intens = 1._X;
                    }
                    else if(during_downramp)
                    {
                        intens = gauss(runTime, 1._X, 0., 1.);
                    }
                    else if(after_endlaser)
                    {
                        intens = 0.;
                    }

                    // add gaussian and sinusoidal pulses

                    for(uint32_t m = 0; m < Unitless::PULSESNR; m++)
                    {
                        float_X intratio = typename Unitless::PULSES_INT_t{}[m];
                        float_X time = typename Unitless::PULSES_TIME_t{}[m] * 1.e-15 / UNIT_TIME;
                        float_X lenratio = typename Unitless::PULSES_LEN_t{}[m];
                        intens += gauss(runTime, intratio, time, lenratio);
                    }

                    for(uint32_t m = 0; m < Unitless::SINPNR; m++)
                    {
                        float_X intratio = typename Unitless::SINPULSES_INT_t{}[m];
                        float_X time1 = typename Unitless::SINPULSES_T1_t{}[m] * 1.e-15 / UNIT_TIME;
                        float_X time2 = typename Unitless::SINPULSES_T2_t{}[m] * 1.e-15 / UNIT_TIME;
                        intens += sinpulse(runTime, intratio, time1, time2);
                    }

                    return Unitless::AMPLITUDE * math::sqrt(intens);
                }

                /** constructor
                 *
                 * @param currentStep current simulation time step
                 */
                HINLINE DivPulses(uint32_t currentStep)
                {
                    // get data
                    DataConnector& dc = Environment<>::get().DataConnector();
                    dataBoxE = dc.get<FieldE>(FieldE::getName(), true)->getDeviceDataBox();

                    // get meta data for offsets
                    SubGrid<simDim> const& subGrid = Environment<simDim>::get().SubGrid();
                    // const DataSpace< simDim > totalCellOffset( subGrid.getGlobalDomain().offset );
                    DataSpace<simDim> const globalCellOffset(subGrid.getLocalDomain().offset);
                    DataSpace<simDim> const halfSimSize(subGrid.getGlobalDomain().size / 2);

                    // transform coordinate system to center of global simulation as origin [cells]
                    offsetToTotalDomain = /* totalCellOffset + */ globalCellOffset - halfSimSize;

                    // @todo reset origin of direction of moving window
                    // offsetToTotalDomain.y() = 0

                    elong = float3_X::create(0.0);

                    /* initialize the laser not in the first cell is equal to a negative shift
                     * in time
                     */
                    const float_64 runTime
                        = Unitless::time_start_init - Unitless::laserTimeShift + DELTA_T * currentStep;

                    phase = float_X(Unitless::w * runTime) + Unitless::LASER_PHASE;

                    float_X const envelope = get_envelope(runTime);

                    if(Unitless::Polarisation == Unitless::LINEAR_X)
                    {
                        elong.x() = envelope * math::sin(phase);
                    }
                    else if(Unitless::Polarisation == Unitless::LINEAR_Z)
                    {
                        elong.z() = envelope * math::sin(phase);
                    }
                    else if(Unitless::Polarisation == Unitless::CIRCULAR)
                    {
                        elong.x() = envelope / math::sqrt(2.0_X) * math::sin(phase);
                        elong.z() = envelope / math::sqrt(2.0_X) * math::cos(phase);
                    }
                }

                /** create device manipulator functor
                 *
                 * @tparam T_WorkerCfg pmacc::mappings::threads::WorkerCfg, configuration of the worker
                 * @tparam T_Acc alpaka accelerator type
                 *
                 * @param alpaka accelerator
                 * @param localSupercellOffset (in supercells, without guards) to the
                 *        origin of the local domain
                 * @param configuration of the worker
                 */
                template<typename T_WorkerCfg, typename T_Acc>
                HDINLINE acc::DivPulses<Unitless> operator()(
                    T_Acc const&,
                    DataSpace<simDim> const& localSupercellOffset,
                    T_WorkerCfg const&) const
                {
                    auto const superCellToLocalOriginCellOffset = localSupercellOffset * SuperCellSize::toRT();
                    return acc::DivPulses<Unitless>(
                        dataBoxE,
                        superCellToLocalOriginCellOffset,
                        offsetToTotalDomain,
                        elong);
                }

                //! get the name of the laser profile
                static HINLINE std::string getName()
                {
                    return "DivPulses";
                }
            };

        } // namespace laserProfiles
    } // namespace fields
} // namespace picongpu
