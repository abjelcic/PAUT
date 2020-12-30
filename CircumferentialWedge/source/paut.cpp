#include <iostream>
#include <cmath>
#include <string>
#include "paut.hh"




Wedge::Wedge(   double          ProbeSeparation                                 ,
                double          PrimaryAxisPitch                                ,
                double          PrimaryAxisSize                                 ,
                double          SecondaryAxisSize                               ,
                unsigned short  NoElementsOnPrimaryAxis                         ,
                double          HeightAtTheMiddleOfTheFirstElement              ,
                double          UltrasonicSpeed                                 ,
                double          PrimaryAxisOffsetOfTheMiddleOfTheFirstElement   ,
                double          SecondaryAxisOffsetOfTheMiddleOfTheFirstElement ,
                double          WedgeWidth                                      ,
                double          RoofAngle                                       ,
                double          WedgeLength                                     ,
                double          WedgeAngle                                      ,
                double          WedgeRadius                                     
             )
{
    constexpr double pi = 4.0 * atan(1.0);
    
    m_N     = NoElementsOnPrimaryAxis;
    m_e     = PrimaryAxisPitch - PrimaryAxisSize;
    m_w     = SecondaryAxisSize;
    m_pitch = PrimaryAxisPitch;
    m_a     = SecondaryAxisOffsetOfTheMiddleOfTheFirstElement - m_w/2.0;
    m_b     = PrimaryAxisOffsetOfTheMiddleOfTheFirstElement - (m_pitch-m_e)/2.0;
    m_omega = WedgeAngle/180.0 * pi;
    m_roof  = RoofAngle /180.0 * pi;
    m_L     = WedgeLength;
    m_W     = WedgeWidth;
    m_R     = WedgeRadius;
    m_c     = UltrasonicSpeed;

    m_d = -2.0*m_W + ProbeSeparation + 2.0*std::cos(m_roof)*( m_a + m_w/2.0 );
    m_H = + HeightAtTheMiddleOfTheFirstElement
        - m_R * ( 1.0 - std::cos( std::asin((m_W+m_d/2.0)/m_R) ) )
        - std::sin(m_roof)  * SecondaryAxisOffsetOfTheMiddleOfTheFirstElement
        - std::sin(m_omega) * PrimaryAxisOffsetOfTheMiddleOfTheFirstElement;

    if( m_e<0.0 || m_a<0.0 || m_b<0.0 || m_d<0.0 || m_H<0.0 )
        throw std::logic_error("Wedge parameters non-physical.");

    return;
}




Specimen::Specimen( double UltrasonicSpeed
                   )
{
    m_c = UltrasonicSpeed;

    return;
}




LawFileParameters::LawFileParameters(   std::string     Version,
                                        unsigned long   Frequency,
                                        unsigned long   Cycles,
                                        unsigned long   SumGain,
                                        unsigned long   Mode,
                                        unsigned long   Filter,
                                        unsigned long   T_First,
                                        unsigned long   R_First,
                                        unsigned long   Scan_Offset,
                                        unsigned long   Index_Offset,
                                        unsigned long   FL_Gain,
                                        unsigned long   Amplitude,
                                        unsigned long   P_width
                                     )
{
    m_Version       = Version;
    m_Frequency     = Frequency;
    m_Cycles        = Cycles;
    m_SumGain       = SumGain;
    m_Mode          = Mode;
    m_Filter        = Filter;
    m_T_First       = T_First;
    m_R_First       = R_First;
    m_Scan_Offset   = Scan_Offset;
    m_Index_Offset  = Index_Offset;
    m_FL_Gain       = FL_Gain;
    m_Amplitude     = Amplitude;
    m_P_width       = P_width;

    return;
}




std::map< SectorScanParameters::FocusingType , std::string > SectorScanParameters::FocusingTypeNames = 
{ 
    { SectorScanParameters::FocusingType::HalfPath     , "HalfPath"     } ,
    { SectorScanParameters::FocusingType::TrueDepth    , "TrueDepth"    } ,
    { SectorScanParameters::FocusingType::AutoFocusing , "AutoFocusing" }
};

SectorScanParameters::SectorScanParameters( SectorScanParameters::FocusingType  FocType,
                                            double                              RefractedAngleStart,
                                            double                              RefractedAngleEnd,
                                            double                              RefractedAngleResolution
                                           )
{
    constexpr double pi = 4.0 * atan(1.0);

    m_FocusingType              = FocType;
    m_RefractedAngleStart       = RefractedAngleStart       / 180.0 * pi;
    m_RefractedAngleEnd         = RefractedAngleEnd         / 180.0 * pi;
    m_RefractedAngleResolution  = RefractedAngleResolution  / 180.0 * pi;


    if( FocType != SectorScanParameters::FocusingType::AutoFocusing )
        throw std::runtime_error("Only supported type of focusing: Auto focus.");

    if( std::abs(RefractedAngleStart) > 90.0    ||
        std::abs(RefractedAngleEnd) > 90.0    ||
        RefractedAngleStart > RefractedAngleEnd ||
        RefractedAngleResolution < 0.0
        )
        throw std::logic_error("Invalid refracted angle boundaries.");

    return;
}




std::ostream & operator << ( std::ostream & out , Wedge const & wedge )
{
    constexpr double pi = 4.0 * atan(1.0);

    out << "Wedge " << std::endl;
    out << "["      << std::endl;
    out << "\t\t"   << "m_N     : " << "\t" << wedge.m_N              << " []" << std::endl;
    out << "\t\t"   << "m_e     : " << "\t" << wedge.m_e              << " [mm]" << std::endl;
    out << "\t\t"   << "m_w     : " << "\t" << wedge.m_w              << " [mm]" << std::endl;
    out << "\t\t"   << "m_pitch : " << "\t" << wedge.m_pitch          << " [mm]" << std::endl;
    out << "\t\t"   << "m_a     : " << "\t" << wedge.m_a              << " [mm]" << std::endl;
    out << "\t\t"   << "m_b     : " << "\t" << wedge.m_b              << " [mm]" << std::endl;
    out << "\t\t"   << "m_omega : " << "\t" << wedge.m_omega/pi*180.0 << " [deg]" << std::endl;
    out << "\t\t"   << "m_roof  : " << "\t" << wedge.m_roof/pi*180.0  << " [deg]" << std::endl;
    out << "\t\t"   << "m_L     : " << "\t" << wedge.m_L              << " [mm]" << std::endl;
    out << "\t\t"   << "m_W     : " << "\t" << wedge.m_W              << " [mm]" << std::endl;
    out << "\t\t"   << "m_R     : " << "\t" << wedge.m_R              << " [mm]" << std::endl;
    out << "\t\t"   << "m_c     : " << "\t" << wedge.m_c              << " [m/s]" << std::endl;
    out << "\t\t"   << "m_d     : " << "\t" << wedge.m_d              << " [mm]" << std::endl;
    out << "\t\t"   << "m_H     : " << "\t" << wedge.m_H              << " [mm]" << std::endl;
    out << "]"      << std::endl;

    return out;
}

std::ostream & operator << ( std::ostream & out , Specimen const & specimen )
{
    out << "Specimen " << std::endl;
    out << "["         << std::endl;
    out << "\t\t"      << "m_c : " << "\t" << specimen.m_c << " [m/s]" << std::endl;
    out << "]"         << std::endl;

    return out;
}

std::ostream & operator << ( std::ostream & out , LawFileParameters const & lawFileParameters )
{
    out << "LawFileParameters " << std::endl;
    out << "["                  << std::endl;
    out << "\t\t"               << "m_Version       : " << "\t" << lawFileParameters.m_Version      << ""       << std::endl;
    out << "\t\t"               << "m_Frequency     : " << "\t" << lawFileParameters.m_Frequency    << " [kHz]" << std::endl;
    out << "\t\t"               << "m_Cycles        : " << "\t" << lawFileParameters.m_Cycles       << " []"    << std::endl;
    out << "\t\t"               << "m_SumGain       : " << "\t" << lawFileParameters.m_SumGain      << " [dB]"  << std::endl;
    out << "\t\t"               << "m_Mode          : " << "\t" << lawFileParameters.m_Mode         << " []"    << std::endl;
    out << "\t\t"               << "m_Filter        : " << "\t" << lawFileParameters.m_Filter       << " []"    << std::endl;
    out << "\t\t"               << "m_T_First       : " << "\t" << lawFileParameters.m_T_First      << " []"    << std::endl;
    out << "\t\t"               << "m_R_First       : " << "\t" << lawFileParameters.m_R_First      << " []"    << std::endl;
    out << "\t\t"               << "m_Scan_Offset   : " << "\t" << lawFileParameters.m_Scan_Offset  << " [um]"  << std::endl;
    out << "\t\t"               << "m_Index_Offset  : " << "\t" << lawFileParameters.m_Index_Offset << " [um]"  << std::endl;
    out << "\t\t"               << "m_FL_Gain       : " << "\t" << lawFileParameters.m_FL_Gain      << " [dB]"  << std::endl;
    out << "\t\t"               << "m_Amplitude     : " << "\t" << lawFileParameters.m_Amplitude    << " [V]"   << std::endl;
    out << "\t\t"               << "m_P_width       : " << "\t" << lawFileParameters.m_P_width      << " [ns]"  << std::endl;
    out << "]"                  << std::endl;

    return out;
}

std::ostream & operator << ( std::ostream & out , SectorScanParameters const & sectorScanParameters )
{
    constexpr double pi = 4.0 * atan(1.0);

    out << "SectorScanParameters " << std::endl;
    out << "["                     << std::endl;
    out << "\t\t"                  << "m_FocusingType             : " << "\t" << SectorScanParameters::FocusingTypeNames[sectorScanParameters.m_FocusingType] << ""       << std::endl;
    out << "\t\t"                  << "m_RefractedAngleStart      : " << "\t" << sectorScanParameters.m_RefractedAngleStart/pi*180.0                          << " [deg]" << std::endl;
    out << "\t\t"                  << "m_RefractedAngleEnd        : " << "\t" << sectorScanParameters.m_RefractedAngleEnd/pi*180.0                            << " [deg]" << std::endl;
    out << "\t\t"                  << "m_RefractedAngleResolution : " << "\t" << sectorScanParameters.m_RefractedAngleResolution/pi*180.0                     << " [deg]" << std::endl;
    out << "]"                     << std::endl;

    return out;
}

