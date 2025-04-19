#include <iostream>
#include "rsexh.hpp"


int main( int argc, char* argv[] )
{
   rsexh::RsExh code;

   std::cout << "N: " << code.N << '\n';

   std::vector<int> a(code.K, 1);
   auto s = rsexh::Encode(a, code.mGf);
   auto c = rsexh::CalculateSyndrome(s, code.R, code.mGf);

   auto a_dec = rsexh::Decode(s, code.R, code.mGf);

   rsexh::show_vector(a, "RS input:");
   rsexh::show_vector(s, "RS output:");
   rsexh::show_vector(c, "cyndrome:");
   rsexh::show_vector(a_dec, "RS decoded:");
   return 0;
}
