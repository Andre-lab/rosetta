// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#include <utility/io/ozstream.hh>
#include <iostream>
#include <fstream>

int
main( int argc, char * argv [] )
{

	utility::io::ozstream(output);
	output.open("bus_error3");
	output << "TEST1\n"  << std::endl;
	output.close();


	std::ofstream output2 ("bus_error.std3");
	output2 << "TEST2\n";
	output2.close();

}
