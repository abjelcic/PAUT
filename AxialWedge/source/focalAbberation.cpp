#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <algorithm>
#include "paut.hh"
#include "utility.hh"
#include "root_rc.hh"
#include "debug.hh"



namespace
{
    std::tuple<double,double,double> generateInitialGuess(  Wedge        const & wedge                ,
                                                            Specimen     const & specimen             ,
                                                            double               RefractedAngle       ,
                                                            double               FocalDepth           ,
                                                            std::size_t          HuygensPrimaryAxis   ,
                                                            std::size_t          HuygensSecondaryAxis 
                                                          )
    {
        auto const SourcePoints = sourcePoints( wedge , HuygensPrimaryAxis , HuygensSecondaryAxis );
        
        Vec<double,3> Savg ( { 0 , 0 , 0 } );
        for( auto const & Si : SourcePoints )
            for( auto const & S : Si )
                Savg += S;
        Savg *= ( 1.0 / double(HuygensPrimaryAxis*HuygensSecondaryAxis*wedge.m_N) );



        double alpha = std::asin( std::sin(RefractedAngle) * wedge.m_c/specimen.m_c );

        
        double FocusX = 0.0;
        double FocusY = + FocalDepth*std::tan(RefractedAngle) + Savg(2) + Savg(3)*std::tan(alpha);
        double FocusZ = - FocalDepth;


        return { FocusX , FocusY , FocusZ };
    }


    double median( std::vector< std::vector<double> > const & ReachingTimes ,
                   std::vector<double>                const & TimeDelays
                 )
    {
        if( ReachingTimes.size() != TimeDelays.size() )
            throw std::logic_error( "Error in focalAbberaton: sizes not compatible." );
        
        std::vector<double> T;
        
        for( unsigned i = 0 ; i < TimeDelays.size() ; ++i )
        {
            double td = TimeDelays[i];
            
            for( double t : ReachingTimes[i] )
                T.push_back( t + td );
        }

        std::sort( T.begin() , T.end() );

        std::size_t n = T.size();
        if( n%2 != 0 )
            return T[n/2];
        else
            return 0.5*( T[(n-1)/2] + T[n/2] );
    }


}








std::tuple< double , Vec<double,3> , std::vector<double> > focalAbberation( Wedge       const & wedge                ,
                                                                            Specimen    const & specimen             ,
                                                                            double              RefractedAngle       ,
                                                                            double              FocalDepth           ,
                                                                            std::size_t         HuygensPrimaryAxis   ,
                                                                            std::size_t         HuygensSecondaryAxis 
                                                                          )
{

    constexpr double pi = 4.0 * std::atan(1.0);

    
    // Initial guess for root-finder
    auto [ FocusX , FocusY , FocusZ ] = generateInitialGuess( wedge               , specimen               ,
                                                              RefractedAngle      , FocalDepth             ,
                                                              HuygensPrimaryAxis  , HuygensSecondaryAxis
                                                             );

    // Root-finding step - search for Focus
    double ferr = -1.0;
    double xerr = -1.0;
    double q[9] = {0,0,0,0,0,0,0,0,0};

    constexpr double      RefAngleTol = 1.e-4; //[deg]
    constexpr double      FocusYTol   = 1.e-4; //[mm]
    constexpr std::size_t maxit       = 100;
    
    unsigned it = 0;
    do
    {
        if( it++ == maxit )
            throw std::runtime_error("Root-finder reached maximum number of iterations.");


        double R = refractedAngle( wedge , specimen , Vec<double,3>({FocusX,FocusY,FocusZ}) , HuygensPrimaryAxis , HuygensSecondaryAxis );


        double x    = FocusY;
        double fx   = ( R - RefractedAngle ) / pi * 180.0;

        FocusY = root_rc( x , fx , ferr , xerr , q );

    } while( ferr >= RefAngleTol && xerr >= FocusYTol  );
    
    DEBUG_MSG( "Root-finder converged after " + std::to_string(it) + " iterations." );










    // Calculation of focal abberation
    auto const SourcePoints = sourcePoints( wedge , HuygensPrimaryAxis , HuygensSecondaryAxis );
    
    auto const ExitPoints = exitPoints( wedge               , specimen                              ,
                                        SourcePoints        , Vec<double,3>({FocusX,FocusY,FocusZ}) ,
                                        HuygensPrimaryAxis , HuygensSecondaryAxis
                                       );

    auto const TimeDelays = timeDelays( wedge.m_c    , specimen.m_c ,
                                        SourcePoints , ExitPoints   ,
                                        Vec<double,3>({FocusX,FocusY,FocusZ})
                                       );
    
    auto const ReachingTimes = reachingTimes( wedge.m_c    , specimen.m_c ,
                                              SourcePoints , ExitPoints   ,
                                              Vec<double,3>({FocusX,FocusY,FocusZ})
                                             );
    
    double const DelayedReachingTimeMedian = median( ReachingTimes , TimeDelays );

    double FocalAbberation = 0.;

    for( unsigned i = 0 ; i < wedge.m_N ; ++i )
    {
        for( unsigned j = 0 ; j < HuygensPrimaryAxis*HuygensSecondaryAxis ; ++j )
        {
            auto   const S      = SourcePoints[i][j];
            auto   const E      = ExitPoints[i][j];
            auto   const F      = Vec<double,3>({FocusX,FocusY,FocusZ});
            double const tDelay = TimeDelays[i];

            double tLeft = 0;
            tLeft += DelayedReachingTimeMedian;
            tLeft -= tDelay;
            tLeft -= (S-E).norm()/wedge.m_c;
            if( tLeft < -1.e-4 )
                throw std::runtime_error( "Error in focalAbberation: Ray never left the wedge." );

            auto const ReachingPoint = E + ( (F-E) * (1.0/(F-E).norm()) ) * specimen.m_c * tLeft;

            FocalAbberation += ( ReachingPoint - F ).norm();
        }
    }


    return { FocalAbberation , Vec<double,3>({FocusX,FocusY,FocusZ}) , TimeDelays };

}

