#include <vector>
#include <tuple>
#include <string>
#include "paut.hh"
#include "utility.hh"
#include "local_min_rc.hh"
#include "debug.hh"

std::tuple<double,double,std::vector<double>> optimizeFocus( Wedge       const & wedge                ,
                                                             Specimen    const & specimen             ,
                                                             double              RefractedAngle       ,
                                                             std::size_t         HuygensPrimaryAxis   ,
                                                             std::size_t         HuygensSecondaryAxis
                                                            )
{


    double FocalDepth;


    // Focal-abberation minimization : finding optimum Focal depth
    double a      =   0.0;
    double b      = 500.0;
    int    status =     0;
    double value  = 1e+10;

    constexpr double      FocalDepthTol = 1.e-4; // [mm]
    constexpr std::size_t maxit         = 100;

    unsigned it = 0;
    do
    {
        if( it++ == maxit )
            throw std::runtime_error("Minimization reached maximum number of iterations.");

        FocalDepth = local_min_rc( a , b , status , value );
        

        auto const [ FocalAbberation , Focus , TimeDelays ] = focalAbberation( wedge              , specimen             ,
                                                                               RefractedAngle     , FocalDepth           ,
                                                                               HuygensPrimaryAxis , HuygensSecondaryAxis
                                                                              );           
        Focus;      // unused
        TimeDelays; // unused

        value = FocalAbberation;


    } while( b-a >= FocalDepthTol );
    
    DEBUG_MSG( "Focal-abberation minimization took " + std::to_string(it) + " iterations." );






    auto const [ FocalAbberation , Focus , TimeDelays ] = focalAbberation( wedge              , specimen             ,
                                                                           RefractedAngle     , FocalDepth           ,
                                                                           HuygensPrimaryAxis , HuygensSecondaryAxis
                                                                          );
    FocalAbberation; // unused

    double const GDelay = gDelay( wedge              , specimen             ,
                                  Focus              , TimeDelays           ,
                                  HuygensPrimaryAxis , HuygensSecondaryAxis
                                 );


    
    return { FocalDepth , GDelay , TimeDelays };                                                                 
}