#include <fstream>
#include <iostream>
#include "debug.hh"
#include "paut.hh"

int main( int argc , char * argv[] )
{
    
    try
    {
        std::string dataFile = "./data/data.in";

        std::ifstream dataStream( dataFile );
        if( !dataStream.is_open() )
        {
            throw std::runtime_error( "Opening file: \"" + dataFile + "\" failed." );
        }
        else
        {
            auto [ wedge , specimen , lawFileParameters , sectorScanParameters ] = readData( dataStream );
            dataStream.close();

            DEBUG_MSG( wedge                );
            DEBUG_MSG( specimen             );
            DEBUG_MSG( lawFileParameters    );
            DEBUG_MSG( sectorScanParameters );


            std::string outLawFile = "./out/AxialProbe.law";
            std::ofstream outLawStream( outLawFile , std::ios::out );
	        if( !outLawStream.is_open()  )
                throw std::runtime_error( "Opening file: \"" + outLawFile + "\" failed." );

            writeLawFile( wedge , specimen , lawFileParameters , sectorScanParameters , outLawStream );

            outLawStream.close();


        }
        
    }
    catch( std::exception const & e )
    {
        std::cerr << e.what() << std::endl;
    }
    



    return 0;
}


