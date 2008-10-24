// $Id: streamtest.cc 4696 2006-08-22 10:01:36Z christi $

/*

  Test to check if the standard streams in libdune can be properly
  linked with this program and if they work

*/

#include "config.h"
#include <fstream>

#include <dune/common/stdstreams.hh>

// enums are a nice special case (was a bug)
enum check { VALUE = 5 };

int main () {
  try {
    // let output happen but vanish
    std::ofstream dummy("/dev/null");   
    Dune::derr.attach(dummy);

    Dune::derr.push(true);
    Dune::derr << "Teststring" << std::endl;

    Dune::derr << VALUE << std::endl;
    Dune::dverb << VALUE << std::endl;
    Dune::dvverb << VALUE << std::endl;
    Dune::dinfo << VALUE << std::endl;
    Dune::dwarn << VALUE << std::endl;
    Dune::dgrave << VALUE << std::endl;
    
    // instantiate private stream and connect global stream
    {
      Dune::DebugStream<> mystream(dummy);
      Dune::derr.tie(mystream);
      Dune::derr << "Blah" << std::endl;
      // untie before mystream gets destructed
      Dune::derr.untie();
    }

    Dune::derr << "Still working" << std::endl;    
  } catch (Dune::Exception &e) {
    std::cerr << e << std::endl;
    return 2;
  } catch (...) {
    return 1;
  };

  return 0;
}