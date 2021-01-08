#include <vector>
#include "paut.hh"
#include "utility.hh"


std::vector<std::vector<double>> reachingTimes( double                                          c_wedge      ,
                                                double                                          c_specimen   ,
                                                std::vector<std::vector<Vec<double,3>>> const & SourcePoints ,
                                                std::vector<std::vector<Vec<double,3>>> const & ExitPoints   ,
                                                Vec<double,3>                           const & Focus
                                                )
{

    bool assert = true;

    assert = assert && ( SourcePoints.size() == ExitPoints.size() );
    assert = assert && ( SourcePoints.size() != 0 );
    for(unsigned i=0; i<SourcePoints.size(); ++i)
        assert = assert && ( SourcePoints[i].size() == ExitPoints[i].size() ) && ( SourcePoints[i].size() != 0 );
    
    if( assert == false )
        throw std::logic_error( "Assertion failed in function reachingTime." );


    std::size_t N = SourcePoints.size();
    std::size_t n = SourcePoints[0].size();

    
    std::vector<std::vector<double>> ReachingTimes(N);
    for( unsigned i = 0 ; i < N ; ++i )
    {
        std::vector<double> T(n);

        for( unsigned j = 0 ; j < n ; ++j )
        {
            auto const S = SourcePoints[i][j];
            auto const E = ExitPoints[i][j];
            auto const F = Focus;

            double t = (S-E).norm()/c_wedge + (E-F).norm()/c_specimen;

            T[j] = t;
        }

        ReachingTimes[i] = T;
    }


    return ReachingTimes;
}