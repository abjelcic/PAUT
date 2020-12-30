#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>
#include "paut.hh"
#include "utility.hh"



std::vector<double> timeDelays( double                                          c_wedge      ,
                                double                                          c_specimen   ,
                                std::vector<std::vector<Vec<double,3>>> const & SourcePoints ,
                                std::vector<std::vector<Vec<double,3>>> const & ExitPoints   ,
                                Vec<double,3>                           const & Focus
                               )
{
    constexpr double pi = 4.0 * std::atan(1.0);
    constexpr std::complex<double> ImagUnit ( 0.0 , 1.0 );


    bool assert = true;

    assert = assert && ( SourcePoints.size() == ExitPoints.size() );
    assert = assert && ( SourcePoints.size() != 0 );
    for(unsigned i=0; i<SourcePoints.size(); ++i)
        assert = assert && ( SourcePoints[i].size() == ExitPoints[i].size() ) && ( SourcePoints[i].size() != 0 );
    
    if( assert == false )
        throw std::logic_error( "Assertion failed in function timeDelays." );





    // Finding reaching times 
    auto const ReachingTimes = reachingTimes( c_wedge , c_specimen , SourcePoints , ExitPoints , Focus );   



    // Finding pinvA and bvec
    std::size_t N = SourcePoints.size();
    std::size_t n = SourcePoints[0].size();

    std::vector< std::vector< std::complex<double> > > Q(N);
    for( unsigned i = 0 ; i < N ; ++i )
    {
        std::vector< std::complex<double> > q(N);
        for( unsigned j = 0 ; j < N ; ++j )
            q[j] = std::exp( ImagUnit * 2.0 * pi/double(N) * double((i+1)*j) ) / std::sqrt(N);

        Q[i] = q; 
    }

    std::vector< std::vector< std::complex<double> > > pinvA;
    pinvA.resize(N);
    for(unsigned i=0; i<N; ++i)
        pinvA[i].resize(N);
    
    for( unsigned row = 0 ; row < N ; ++row )
        for( unsigned col = 0 ; col < N ; ++col )
        {
            std::complex<double> z {0};
            for( unsigned k = 0 ; k < N-1 ; ++k )
                z += Q[k][row] * std::conj(Q[k][col]);

            pinvA[row][col] = 1.0/double(n) * z;
        }

    std::vector< std::complex<double> > bvec(N);
    for( unsigned i = 0 ; i < N ; ++i )
    {
        std::complex<double> b {0};

        for( unsigned j = 0 ; j < N ; ++j )
            for( unsigned k = 0 ; k < n ; ++k )
                b += ReachingTimes[j][k] * ( 1.0/double(N) - ( (i==j) ? 1.0 : 0.0 ) );

        bvec[i] = b;
    }



    // Finding time delays
    std::vector<double> TimeDelays(N);
    for( unsigned i = 0 ; i < N ; ++i )
    {
        double td = 0;
        for( unsigned j = 0 ; j < N ; ++j )
            td += std::real( pinvA[i][j] * bvec[j] );

        TimeDelays[i] = td;
    }

    
    double const min = *std::min_element( TimeDelays.begin() , TimeDelays.end() );
    std::for_each( TimeDelays.begin() , TimeDelays.end() , [min](double &x){ x -= min; } );



    return TimeDelays;
}


