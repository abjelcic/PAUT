#include <vector>
#include "paut.hh"
#include "utility.hh"

std::vector< std::vector< Vec<double,3> > > sourcePoints( Wedge       const & wedge                 ,    
                                                          std::size_t         HuygensPrimaryAxis    ,
                                                          std::size_t         HuygensSecondaryAxis
                                                        )

{

    std::vector< std::vector< Vec<double,3> > > Sources;
    Sources.reserve( wedge.m_N );
    for( unsigned i = 1 ; i <= wedge.m_N ; ++i )
    {
        double omega = wedge.m_omega;
        double roof  = wedge.m_roof;
        double W     = wedge.m_W;
        double L     = wedge.m_L;
        double d     = wedge.m_d;
        double H     = wedge.m_H;
        double R     = wedge.m_R;
        double a     = wedge.m_a;
        double b     = wedge.m_b;
        double pitch = wedge.m_pitch;
        double w     = wedge.m_w;
        double e     = wedge.m_e;


        Vec<double,3> const e1 ( {               0.0 ,
                                   + std::cos(omega) ,
                                   + std::sin(omega)
                                  }
                                );
        Vec<double,3> const e2 ( { - std::cos(roof) ,
                                                0.0 ,
                                   + std::sin(roof) ,
                                  }
                                );

        auto const P =   Vec<double,3>( { W + 0.5*d , 0.0 , H + R - std::sqrt( R*R - (L/2)*(L/2) ) } )
                       + ( b+(i-1)*pitch ) * e1
                       +                a  * e2;


        std::vector< Vec<double,3> > Sources_i;
        Sources_i.reserve( HuygensPrimaryAxis * HuygensSecondaryAxis );

        for( unsigned i1 = 1 ; i1 <= HuygensPrimaryAxis ; ++i1 )
        {
            for( unsigned i2 = 1 ; i2 <= HuygensSecondaryAxis ; ++i2 )
            {
                double N1 = double(HuygensPrimaryAxis  );
                double N2 = double(HuygensSecondaryAxis);

                auto const S = P + ( (i1-0.5)*(pitch-e)/N1 )*e1 + ( (i2-0.5)*w/N2 )*e2;

                Sources_i.push_back( S );
            }
        }
        
        Sources.push_back( Sources_i );
    }

    return Sources;
}