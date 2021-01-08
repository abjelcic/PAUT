#include <vector>
#include "paut.hh"
#include "utility.hh"



std::vector<std::vector<Vec<double,3>>> exitPoints( Wedge                                   const & wedge               ,
                                                    Specimen                                const & specimen             ,
                                                    std::vector<std::vector<Vec<double,3>>> const & SourcePoints         ,
                                                    Vec<double,3>                           const & Focus                ,
                                                    std::size_t                                     HuygensPrimaryAxis   ,
                                                    std::size_t                                     HuygensSecondaryAxis
                                                    )
{

    auto CylLineIntersection = [] ( Wedge::WedgeType wedgeType , double R , double L , Vec<double,3> In , Vec<double,3> Out ) -> Vec<double,3>
    {
        #ifdef DEBUG
        if( wedgeType == Wedge::WedgeType::Axial )
        {
            if( In(1)*In(1) + (In(3)-R)*(In(3)-R) > R*R )
                throw std::logic_error("CylLineIntersection: \"In\" vector is in the exterior of the cylinder.");
            if( Out(1)*Out(1) + (Out(3)-R)*(Out(3)-R) < R*R )
                throw std::logic_error("CylLineIntersection: \"Out\" vector is in the interior of the cylinder.");
        }
        else
        {
            if( (In(2)-L/2)*(In(2)-L/2) + (In(3)-R)*(In(3)-R) > R*R )
                throw std::logic_error("CylLineIntersection: \"In\" vector is in the exterior of the cylinder.");
            if( (Out(2)-L/2)*(Out(2)-L/2) + (Out(3)-R)*(Out(3)-R) < R*R )
                throw std::logic_error("CylLineIntersection: \"Out\" vector is in the interior of the cylinder.");            
        }
        #endif

        auto const v = Out - In;
        double a;
        double b;
        double c;

        if( wedgeType == Wedge::WedgeType::Axial )
        {
            a = v(1)*v(1) + v(3)*v(3);
            b = 2*( In(1)*v(1) + In(3)*v(3) - 2*R*v(3) );
            c = In(1)*In(1) + (In(3)-R)*(In(3)-R) - R*R;
        }
        else
        {
            a = v(2)*v(2) + v(3)*v(3);
            b = 2*( (In(2)-L/2)*v(2) + (In(3)-R)*v(3) );
            c = (In(2)-L/2)*(In(2)-L/2) + (In(3)-R)*(In(3)-R) - R*R;            
        }

        double t = ( -b + std::sqrt( b*b - 4*a*c ) ) / ( 2*a );

        return In + t*v;
    };

    

    std::vector< std::vector< Vec<double,3> > > ExitPoints ( SourcePoints.size() );
    for( unsigned i = 0 ; i < SourcePoints.size() ; ++i )
    {
        std::vector< Vec<double,3> > ExitPoints_i ( SourcePoints[i].size() );
        for( unsigned j = 0 ; j < SourcePoints[i].size(); ++j )
        {
            auto const S      = SourcePoints[i][j];
            auto const F      = Focus;
            auto const Eguess = ( j!=0 ) ? ExitPoints_i[j-1] : CylLineIntersection( wedge.m_WedgeType , wedge.m_R , wedge.m_L , S , F );

            auto const E      = solveCylSnell( wedge.m_WedgeType , wedge.m_c , specimen.m_c , wedge.m_R , wedge.m_L , S , F , Eguess );

            ExitPoints_i[j]   = E;
        }

        ExitPoints[i] = ExitPoints_i;
    }


    return ExitPoints;

}


