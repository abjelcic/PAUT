#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include "paut.hh"
#include "utility.hh"


void writeLawFile( Wedge                const & wedge                ,
                   Specimen             const & specimen             ,
                   LawFileParameters    const & lawFileParameters    ,
                   SectorScanParameters const & sectorScanParameters ,
                   std::ofstream              & outLawStream            
                  )
{

    auto RefractedAngles = uniformGrid1D( sectorScanParameters.m_RefractedAngleStart      ,
                                          sectorScanParameters.m_RefractedAngleResolution ,
                                          sectorScanParameters.m_RefractedAngleEnd        ,
                                          1.e-3
                                         );
    
    std::size_t N_laws = RefractedAngles.size();

    outLawStream << "V" + lawFileParameters.m_Version + " " + std::to_string(N_laws) << std::endl;


    
    for( unsigned i = 0 ; i < N_laws ; ++i )
    {
        auto R = RefractedAngles[i];

        switch( sectorScanParameters.m_FocusingType )
        {
            case SectorScanParameters::FocusingType::AutoFocusing :
            {
                auto const [ FocalDepth , GDelay , TimeDelays ] = optimizeFocus( wedge , specimen , R );
                
                constexpr double pi = 4.0 * std::atan(1.0);

                outLawStream <<     wedge.m_N                           << " ";
                outLawStream <<     lawFileParameters.m_Frequency       << " ";
                outLawStream <<     lawFileParameters.m_Cycles          << " ";
                outLawStream <<     lawFileParameters.m_SumGain         << " ";
                outLawStream <<     lawFileParameters.m_Mode            << " ";
                outLawStream <<     lawFileParameters.m_Filter          << " ";
                outLawStream <<     std::round(   R/pi*180.0 * 10.0 )   << " ";
                outLawStream <<     std::round( 0.0/pi*180.0 * 10.0 )   << " ";
                outLawStream <<     lawFileParameters.m_T_First         << " ";
                outLawStream <<     lawFileParameters.m_R_First         << " ";
                outLawStream <<     lawFileParameters.m_Scan_Offset     << " ";
                outLawStream <<     lawFileParameters.m_Index_Offset    << " ";
                outLawStream <<     std::round( GDelay * 1.e+6 )        << " ";
                outLawStream <<     std::round( FocalDepth * 1.e+3 )    << " ";
                outLawStream <<     specimen.m_c                        << " ";
                outLawStream <<     0                                   << std::endl;

                assert( wedge.m_N == TimeDelays.size() );

                for( unsigned j = 0 ; j < TimeDelays.size() ; ++j )
                {
                    outLawStream << j+1                                 << " ";
                    outLawStream << lawFileParameters.m_FL_Gain         << " ";
                    outLawStream << std::round( TimeDelays[j] * 1.e+6 ) << " ";
                    outLawStream << std::round( TimeDelays[j] * 1.e+6 ) << " ";
                    outLawStream << lawFileParameters.m_Amplitude       << " ";
                    outLawStream << lawFileParameters.m_P_width         <<  "";
                    
                    if( !( j==TimeDelays.size()-1 && i==N_laws-1 ) )
                        outLawStream << std::endl;    
                }

                break;
            }
            default :
            {
                throw std::logic_error("Only Auto focus supported as Focusing type.");
                break;
            }
        }





    }

    return;
}
