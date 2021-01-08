#include <cmath>
#include <utility>
#include <string>
#include "paut.hh"
#include "utility.hh"
#include "optmin.hh"
#include "debug.hh"


Vec<double,3> solveCylSnellCircumferential( double                c_wedge    ,
                                            double                c_specimen ,
                                            double                Radius     ,
                                            double                Length     ,
                                            Vec<double,3> const & Source     ,
                                            Vec<double,3> const & Focus      ,
                                            Vec<double,3> const & ExitGuess  
                                           );

Vec<double,3> solveCylSnellAxial( double                c_wedge    ,
                                  double                c_specimen ,
                                  double                Radius     ,
                                  double                Length     ,
                                  Vec<double,3> const & Source     ,
                                  Vec<double,3> const & Focus      ,
                                  Vec<double,3> const & ExitGuess  
                                 );



Vec<double,3> solveCylSnellCircumferential( double                c_wedge    ,
                                            double                c_specimen ,
                                            double                Radius     ,
                                            double                Length     ,
                                            Vec<double,3> const & Source     ,
                                            Vec<double,3> const & Focus      ,
                                            Vec<double,3> const & ExitGuess  
                                           )
{

    auto CartesianToCylindrical = []( Vec<double,3> const & T , double R , double L ) -> std::pair<double,double>
    {
        double x0 = T(1);
        double y0 = T(2);
        double z0 = T(3);

        double x = R - z0;
        double y = y0 - L/2;
        double z = x0;

        return { z , std::atan2(y,x) };
    };

    auto [ x0 , phi0 ] = CartesianToCylindrical( ExitGuess , Radius , Length );




    auto f = [=]( double xmin[] ) -> double
    {
        double x   = xmin[0];
        double phi = xmin[1];
        double R   = Radius;
        double L   = Length;

        double ans = + ( Vec<double,3>( { x , L/2 + R*std::sin(phi) , R*(1-std::cos(phi)) } ) - Source ).norm() / c_wedge
                     + ( Vec<double,3>( { x , L/2 + R*std::sin(phi) , R*(1-std::cos(phi)) } ) - Focus  ).norm() / c_specimen;

        return ans;
    };
    
    constexpr std::size_t n = 2;

    double start[2] = { x0 , phi0 };

    double xmin[2] = {0};
    
    double fmin = -1;

    double reqmin = 1.e-5;

    double step[2] = { 2 , 0.1 };

    int konvgecheck = 25;

    int maxfuncalls = 1000;

    int icount = -1;

    int numrestarts = -1;

    int ifault = -1;

    nelmin( f           , n           ,
            start       , xmin        , &fmin   ,
            reqmin      , step        ,
            konvgecheck , maxfuncalls , &icount , &numrestarts ,
            &ifault
           );
    
    DEBUG_MSG( "Nelder-Mead error indicator: " + std::to_string(ifault) +
               ", number of function calls: "  + std::to_string(icount)   );
    

    double   xTarget = xmin[0];
    double phiTarget = xmin[1];

    return Vec<double,3>( { xTarget                               ,
                            Length/2 + Radius*std::sin(phiTarget) ,
                            Radius*(1-std::cos(phiTarget))
                           }
                         );

}



Vec<double,3> solveCylSnellAxial( double                c_wedge    ,
                                  double                c_specimen ,
                                  double                Radius     ,
                                  double                Length     ,
                                  Vec<double,3> const & Source     ,
                                  Vec<double,3> const & Focus      ,
                                  Vec<double,3> const & ExitGuess  
                                 )
{

    auto CartesianToCylindrical = []( Vec<double,3> const & T , double R ) -> std::pair<double,double>
    {
        double x0 = T(1);
        double y0 = T(2);
        double z0 = T(3);

        double x = R - z0;
        double y = x0;
        double z = y0;

        return { z , std::atan2(y,x) };
    };

    auto [ z0 , phi0 ] = CartesianToCylindrical( ExitGuess , Radius );




    auto f = [=]( double xmin[] ) -> double
    {
        double z   = xmin[0];
        double phi = xmin[1];
        double R   = Radius;

        double ans = + ( Vec<double,3>( { R*std::sin(phi) , z , R*(1-std::cos(phi)) } ) - Source ).norm() / c_wedge
                     + ( Vec<double,3>( { R*std::sin(phi) , z , R*(1-std::cos(phi)) } ) - Focus  ).norm() / c_specimen;

        return ans;
    };
    
    constexpr std::size_t n = 2;

    double start[2] = { z0 , phi0 };

    double xmin[2] = {0};
    
    double fmin = -1;

    double reqmin = 1.e-5;

    double step[2] = { 2 , 0.1 };

    int konvgecheck = 25;

    int maxfuncalls = 1000;

    int icount = -1;

    int numrestarts = -1;

    int ifault = -1;

    nelmin( f           , n           ,
            start       , xmin        , &fmin   ,
            reqmin      , step        ,
            konvgecheck , maxfuncalls , &icount , &numrestarts ,
            &ifault
           );
    
    DEBUG_MSG( "Nelder-Mead error indicator: " + std::to_string(ifault) +
               ", number of function calls: "  + std::to_string(icount)   );
    

    double   zTarget = xmin[0];
    double phiTarget = xmin[1];

    return Vec<double,3>( { Radius*std::sin(phiTarget)     ,
                            zTarget                        ,
                            Radius*(1-std::cos(phiTarget))
                           }
                         );

}



Vec<double,3> solveCylSnell( Wedge::WedgeType      wedgeType  ,
                             double                c_wedge    ,
                             double                c_specimen ,
                             double                Radius     ,
                             double                Length     ,
                             Vec<double,3> const & Source     ,
                             Vec<double,3> const & Focus      ,
                             Vec<double,3> const & ExitGuess  
                            )
{

    switch( wedgeType )
    {
        case Wedge::WedgeType::Axial :
            
            return solveCylSnellAxial( c_wedge , c_specimen , Radius , Length , Source , Focus , ExitGuess );
            break;

        case Wedge::WedgeType::Circumferential :
            
            return solveCylSnellCircumferential( c_wedge , c_specimen , Radius , Length , Source , Focus , ExitGuess );
            break;

        default:
            throw std::logic_error("Non-existing wedge type in solveCylSnell()!");
            break;
    }

}

