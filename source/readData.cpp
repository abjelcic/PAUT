#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <tuple>
#include "paut.hh"

std::tuple<Wedge,Specimen,LawFileParameters,SectorScanParameters> readData( std::ifstream & DataFile )
{

    std::vector<std::string> lines;
    for( std::string line ; std::getline( DataFile , line ) ; )
        if( line.size() > 0)
            lines.emplace_back(line);


    Wedge::WedgeType wedgeType;
    double           ProbeSeparation;
    double           PrimaryAxisPitch;
    double           PrimaryAxisSize;
    double           SecondaryAxisSize;
    unsigned short   NoElementsOnPrimaryAxis;
    double           RefractedAngleStart;
    double           RefractedAngleEnd;
    double           RefractedAngleResolution;
    double           SpecimenUltrasoundSpeed;
    double           HeightAtTheMiddleOfTheFirstElement;
    double           WedgeUltrasonicSpeed;
    double           PrimaryAxisOffsetOfTheMiddleOfTheFirstElement;
    double           SecondaryAxisOffsetOfTheMiddleOfTheFirstElement;
    double           WedgeWidth;
    double           RoofAngle;
    double           WedgeLength;
    double           WedgeAngle;
    double           WedgeRadius;
    std::string      Version;
    unsigned long    Frequency;
    unsigned long    Cycles;
    unsigned long    SumGain;
    unsigned long    Mode;
    unsigned long    Filter;
    unsigned long    T_First;
    unsigned long    R_First;
    unsigned long    Scan_Offset;
    unsigned long    Index_Offset;
    unsigned long    FL_Gain;
    unsigned long    Amplitude;
    unsigned long    P_width;
    
    if      ( lines[0] == std::string("Axial")           )
    {
        wedgeType = Wedge::WedgeType::Axial;
    }
    else if ( lines[0] == std::string("Circumferential") )
    {
        wedgeType = Wedge::WedgeType::Circumferential;
    }
    else
    {
        throw std::logic_error("Wedge type must be \"Axial\" or \"Circumferential\"!");
    }

    unsigned int i = 1;
    std::stringstream( lines[i++] ) >> ProbeSeparation;
    std::stringstream( lines[i++] ) >> PrimaryAxisPitch;
    std::stringstream( lines[i++] ) >> PrimaryAxisSize;
    std::stringstream( lines[i++] ) >> SecondaryAxisSize;
    std::stringstream( lines[i++] ) >> NoElementsOnPrimaryAxis;
    std::stringstream( lines[i++] ) >> RefractedAngleStart;
    std::stringstream( lines[i++] ) >> RefractedAngleEnd;
    std::stringstream( lines[i++] ) >> RefractedAngleResolution;
    std::stringstream( lines[i++] ) >> SpecimenUltrasoundSpeed;
    std::stringstream( lines[i++] ) >> HeightAtTheMiddleOfTheFirstElement;
    std::stringstream( lines[i++] ) >> WedgeUltrasonicSpeed;
    std::stringstream( lines[i++] ) >> PrimaryAxisOffsetOfTheMiddleOfTheFirstElement;
    std::stringstream( lines[i++] ) >> SecondaryAxisOffsetOfTheMiddleOfTheFirstElement;
    std::stringstream( lines[i++] ) >> WedgeWidth;
    std::stringstream( lines[i++] ) >> RoofAngle;
    std::stringstream( lines[i++] ) >> WedgeLength;
    std::stringstream( lines[i++] ) >> WedgeAngle;
    std::stringstream( lines[i++] ) >> WedgeRadius;
    std::stringstream( lines[i++] ) >> Version;
    std::stringstream( lines[i++] ) >> Frequency;
    std::stringstream( lines[i++] ) >> Cycles;
    std::stringstream( lines[i++] ) >> SumGain;
    std::stringstream( lines[i++] ) >> Mode;
    std::stringstream( lines[i++] ) >> Filter;
    std::stringstream( lines[i++] ) >> T_First;
    std::stringstream( lines[i++] ) >> R_First;
    std::stringstream( lines[i++] ) >> Scan_Offset;
    std::stringstream( lines[i++] ) >> Index_Offset;
    std::stringstream( lines[i++] ) >> FL_Gain;
    std::stringstream( lines[i++] ) >> Amplitude;
    std::stringstream( lines[i++] ) >> P_width;
    




    Wedge wedge( wedgeType,
                 ProbeSeparation,
                 PrimaryAxisPitch,
                 PrimaryAxisSize,
                 SecondaryAxisSize,
                 NoElementsOnPrimaryAxis,
                 HeightAtTheMiddleOfTheFirstElement,
                 WedgeUltrasonicSpeed,
                 PrimaryAxisOffsetOfTheMiddleOfTheFirstElement,
                 SecondaryAxisOffsetOfTheMiddleOfTheFirstElement,
                 WedgeWidth,
                 RoofAngle,
                 WedgeLength,
                 WedgeAngle,
                 WedgeRadius
                );

    Specimen specimen( SpecimenUltrasoundSpeed );

    SectorScanParameters sectorScanParameters( SectorScanParameters::FocusingType::AutoFocusing,
                                               RefractedAngleStart,   
                                               RefractedAngleEnd,       
                                               RefractedAngleResolution 
                                              );

    LawFileParameters lawFileParameters( Version,
                                         Frequency,
                                         Cycles,
                                         SumGain,
                                         Mode,
                                         Filter,
                                         T_First,
                                         R_First,
                                         Scan_Offset,
                                         Index_Offset,
                                         FL_Gain,
                                         Amplitude,
                                         P_width
                                        );


    return { wedge , specimen , lawFileParameters , sectorScanParameters };

}


