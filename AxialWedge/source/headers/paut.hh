#ifndef PAUT_H
#define PAUT_H

#include <iostream>
#include <string>
#include <fstream>
#include <tuple>
#include <map>
#include <vector>
#include "utility.hh"

struct Wedge
{
    unsigned short  m_N;        //[]
    double          m_e;        //[mm]
    double          m_w;        //[mm]
    double          m_pitch;    //[mm]
    double          m_a;        //[mm]
    double          m_b;        //[mm]
    double          m_omega;    //[rad]
    double          m_roof;     //[rad]
    double          m_L;        //[mm]
    double          m_W;        //[mm]
    double          m_R;        //[mm]
    double          m_c;        //[m/s]
    double          m_d;        //[mm]
    double          m_H;        //[mm]

    Wedge(  double          ProbeSeparation                                 , //[mm]
            double          PrimaryAxisPitch                                , //[mm]
            double          PrimaryAxisSize                                 , //[mm]
            double          SecondaryAxisSize                               , //[mm]
            unsigned short  NoElementsOnPrimaryAxis                         , //[]
            double          HeightAtTheMiddleOfTheFirstElement              , //[mm]
            double          UltrasonicSpeed                                 , //[m/s]
            double          PrimaryAxisOffsetOfTheMiddleOfTheFirstElement   , //[mm]
            double          SecondaryAxisOffsetOfTheMiddleOfTheFirstElement , //[mm]
            double          WedgeWidth                                      , //[mm]
            double          RoofAngle                                       , //[deg]
            double          WedgeLength                                     , //[mm]
            double          WedgeAngle                                      , //[deg]
            double          WedgeRadius                                       //[mm]
          );
};

struct Specimen
{
    double m_c; // [m/s]

    Specimen( double UltrasonicSpeed //[m/s]
             );
};

struct LawFileParameters
{
    std::string     m_Version;
    unsigned long   m_Frequency;
    unsigned long   m_Cycles;
    unsigned long   m_SumGain;
    unsigned long   m_Mode;
    unsigned long   m_Filter;
    unsigned long   m_T_First;
    unsigned long   m_R_First;
    unsigned long   m_Scan_Offset;
    unsigned long   m_Index_Offset;
    unsigned long   m_FL_Gain;
    unsigned long   m_Amplitude;
    unsigned long   m_P_width;

    LawFileParameters(  std::string     Version,
                        unsigned long   Frequency,
                        unsigned long   Cycles,
                        unsigned long   SumGain,
                        unsigned long   Mode,
                        unsigned long   Filter,
                        unsigned long   T_First,
                        unsigned long   R_First,
                        unsigned long   Scan_Offset,
                        unsigned long   Index_Offset,
                        unsigned long   FL_Gain,
                        unsigned long   Amplitude,
                        unsigned long   P_width
                     );

};

struct SectorScanParameters
{
    enum class FocusingType
    {
        HalfPath,
        TrueDepth,
        AutoFocusing
    };

    static std::map< FocusingType , std::string > FocusingTypeNames;


    FocusingType    m_FocusingType;
    double          m_RefractedAngleStart;      // [rad]
    double          m_RefractedAngleEnd;        // [rad]
    double          m_RefractedAngleResolution; // [rad]

    SectorScanParameters(   SectorScanParameters::FocusingType  FocType,
                            double                              RefractedAngleStart,     // [deg]
                            double                              RefractedAngleEnd,       // [deg]
                            double                              RefractedAngleResolution // [deg]
                        );

};

std::ostream & operator << ( std::ostream & out , Wedge                 const & wedge                );
std::ostream & operator << ( std::ostream & out , Specimen              const & specimen             );
std::ostream & operator << ( std::ostream & out , LawFileParameters     const & lawFileParameters    );
std::ostream & operator << ( std::ostream & out , SectorScanParameters  const & sectorScanParameters );

std::tuple<Wedge,Specimen,LawFileParameters,SectorScanParameters> readData( std::ifstream & DataFile );

void writeLawFile( Wedge                const & wedge                ,
                   Specimen             const & specimen             ,
                   LawFileParameters    const & lawFileParameters    ,
                   SectorScanParameters const & sectorScanParameters ,
                   std::ofstream              & outLawStream            
                  );

Vec<double,3> solveCylSnell( double                c_wedge    ,
                             double                c_specimen ,
                             double                Radius     ,
                             Vec<double,3> const & Source     ,
                             Vec<double,3> const & Focus      ,
                             Vec<double,3> const & ExitGuess  
                            );

double refractedAngle( Wedge         const & wedge                     ,
                       Specimen      const & specimen                  ,
                       Vec<double,3> const & Focus                     ,
                       std::size_t           HuygensPrimaryAxis   = 1  ,
                       std::size_t           HuygensSecondaryAxis = 5
                      );

std::vector< std::vector< Vec<double,3> > > sourcePoints( Wedge       const & wedge                    ,    
                                                          std::size_t         HuygensPrimaryAxis   = 1 ,
                                                          std::size_t         HuygensSecondaryAxis = 5
                                                        );

std::vector<std::vector<Vec<double,3>>> exitPoints( Wedge                                   const & wedge                     ,
                                                    Specimen                                const & specimen                  ,
                                                    std::vector<std::vector<Vec<double,3>>> const & SourcePoints              ,
                                                    Vec<double,3>                           const & Focus                     ,
                                                    std::size_t                                     HuygensPrimaryAxis   = 1  ,
                                                    std::size_t                                     HuygensSecondaryAxis = 5
                                                    );

// returns: { Focal abberation , Focus , Time delays }
std::tuple< double , Vec<double,3> , std::vector<double> > focalAbberation( Wedge       const & wedge                    ,
                                                                            Specimen    const & specimen                 ,
                                                                            double              RefractedAngle           ,
                                                                            double              FocalDepth               ,
                                                                            std::size_t         HuygensPrimaryAxis   = 1 ,
                                                                            std::size_t         HuygensSecondaryAxis = 5
                                                                          );


std::vector<double> timeDelays( double                                          c_wedge      ,
                                double                                          c_specimen   ,
                                std::vector<std::vector<Vec<double,3>>> const & SourcePoints ,
                                std::vector<std::vector<Vec<double,3>>> const & ExitPoints   ,
                                Vec<double,3>                           const & Focus
                               );

std::vector<std::vector<double>> reachingTimes( double                                          c_wedge      ,
                                                double                                          c_specimen   ,
                                                std::vector<std::vector<Vec<double,3>>> const & SourcePoints ,
                                                std::vector<std::vector<Vec<double,3>>> const & ExitPoints   ,
                                                Vec<double,3>                           const & Focus
                                                );       

// returns: { Focal depth , Global delay , Time delays }
std::tuple<double,double,std::vector<double>> optimizeFocus( Wedge       const & wedge                    ,
                                                             Specimen    const & specimen                 ,
                                                             double              RefractedAngle           ,
                                                             std::size_t         HuygensPrimaryAxis   = 1 ,
                                                             std::size_t         HuygensSecondaryAxis = 5
                                                            );

double gDelay( Wedge               const & wedge                    ,
               Specimen            const & specimen                 ,
               Vec<double,3>       const & Focus                    ,
               std::vector<double> const & TimeDelays               ,
               std::size_t                 HuygensPrimaryAxis   = 1 ,
               std::size_t                 HuygensSecondaryAxis = 5
              );


#endif // PAUT_H


