#include <cmath>
#include <vector>
#include <cmath>
#include "paut.hh"
#include "utility.hh"



double refractedAngle( Wedge         const & wedge                ,
                       Specimen      const & specimen             ,
                       Vec<double,3> const & Focus                ,
                       std::size_t           HuygensPrimaryAxis   ,
                       std::size_t           HuygensSecondaryAxis
                      )
{

    auto const Sources = sourcePoints( wedge , HuygensPrimaryAxis , HuygensSecondaryAxis );

    auto const ExitPoints = exitPoints( wedge , specimen , Sources , Focus , HuygensPrimaryAxis , HuygensSecondaryAxis );

    Vec<double,3> V ( { 0 , 0 , 0 } );
    for( auto const & ExitPoints_i : ExitPoints )
        for( auto const & E : ExitPoints_i )
            V += ( Focus - E );
    V(1) = 0.0;
    V *= (1/V.norm()); 

    double RefractedAngle = std::acos( Vec<double,3>::dot( V , Vec<double,3>({0,0,-1}) ) );
    RefractedAngle = ( V(2) < 0 ) ? -RefractedAngle : RefractedAngle;

    return RefractedAngle;

}
