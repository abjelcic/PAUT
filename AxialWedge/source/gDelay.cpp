#include <vector>
#include <numeric>
#include "paut.hh"



namespace
{
    double average( std::vector< std::vector<double> > v )
    {
        double      ans = 0.;
        std::size_t n   = 0;

        for( auto const & vi : v )
            for( auto const & x : vi )
            {
                n++;
                ans += x;
            }
        
        return ans/double(n);
    }
}

double gDelay( Wedge               const & wedge                ,
               Specimen            const & specimen             ,
               Vec<double,3>       const & Focus                ,
               std::vector<double> const & TimeDelays           ,
               std::size_t                 HuygensPrimaryAxis   ,
               std::size_t                 HuygensSecondaryAxis
              )
{

    auto const SourcePoints = sourcePoints( wedge , HuygensPrimaryAxis , HuygensSecondaryAxis );
    
    auto const ExitPoints = exitPoints( wedge               , specimen              ,
                                        SourcePoints        , Focus                 ,
                                        HuygensPrimaryAxis  , HuygensSecondaryAxis
                                       );
    
    double GDelay = 0.0;
    std::size_t n = 0;
    for( unsigned int i = 0 ; i < wedge.m_N ; ++i )
    {
        for( unsigned j = 0 ; j < HuygensPrimaryAxis*HuygensSecondaryAxis ; ++j )
        {
            n++;

            auto const S = SourcePoints[i][j];
            auto const E = ExitPoints[i][j];

            GDelay += 2.0 * ( (S-E).norm()/wedge.m_c + TimeDelays[i] ); 
        }
    }

    GDelay /= double(n);

    return GDelay;
}
 
