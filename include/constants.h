#include <boost/math/constants/constants.hpp>
#include <cmath>

namespace constants {

    /* UNITS USED:
    distance nm
    time     ps
    energy   eV
    */
const double pi = boost::math::constants::pi<double>();

// The fine structure constant (dimensionless)
const double alpha = 0.0072973525693;
// hbar eV*ps
const double hbar = 6.5821192815E-4;
const double hbar2 = std::pow(6.5821192815E-4,2);
const double hbar3 = std::pow(6.5821192815E-4,3);
// speed of light nm / ps
const double c = 299792.458;
const double c2 = std::pow(299792.458, 2); 

/***************************/
/*       CONVERSIONS       */
/***************************/
const double bohr2nm = 0.052917721092;
const double nm2bohr = 18.897259886;
const double ang2bohr = 1.8897259886;
const double bohr2ang = 1.0 / 1.8897259886;
const double nm2ang = 10.0;
const double ang2nm = 0.1;
const double hrt2ev = 27.21138602;
const double ev2hrt = 1.0 / 27.21138602;

} // namespace constants