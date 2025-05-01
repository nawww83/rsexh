#include <iostream>
#include <random>
#include "rsexh.hpp"

static auto const seed = std::random_device{}();

/***
 * Генератор случайных чисел.
 */
auto roll_uint = [urbg = std::mt19937{seed},
                distr = std::uniform_int_distribution<uint>{}]() mutable -> uint {
   return distr(urbg);
};

auto roll_error = [urbg = std::mt19937{seed},
                distr = std::uniform_int_distribution<uint>{}](double error_probability) mutable -> bool {
   const uint threshold = error_probability * uint(-1);
   return error_probability == 0 ? false : distr(urbg) <= threshold;
};

void test_rs() {
   rsexh::RsExh code;

   std::cout << "N: " << code.N << '\n';

   std::vector<int> a(code.K, 14);
   auto s = rsexh::Encode(a, code.mGf);
   auto c = rsexh::CalculateSyndrome(s, code.R, code.mGf);

   auto a_dec = rsexh::Decode(s, code.R, code.mGf);

   rsexh::show_vector(a, "RS input:");
   rsexh::show_vector(s, "RS output:");
   rsexh::show_vector(c, "cyndrome:");
   rsexh::show_vector(a_dec, "RS decoded:");
}

void measure_ber() {
   rsexh::RsExh code;
   std::cout << "N: " << code.N << '\n';
   hamming::CodeWord<int, code.M2> a(code.mHammingCode.K);
   hamming::CodeWord<int, code.M2> a_received(code.mHammingCode.N);
   // Source
   for (auto& el : a) {
      el.mStatus = hamming::SymbolStatus::Normal;
      for (auto& symbol : el.mSymbol)
         symbol = roll_uint() & 15; // Полубайты.
   }
   // hamming::show_codeword(a, code, "Input a: ");
   // Hamming encode
   std::cout << "Hamming encode\n";
   auto s_h = code.mHammingCode.Encode(a);
   std::vector<int> a_rs(code.K);
   std::vector<std::vector<int>> v;
   // RS encode + CRC
   std::cout << "RS encode\n";
   v.clear();
   for (const auto& el : s_h) {
      assert(code.K == el.mSymbol.size());
      for (int i = 0; i < code.K; ++i)
         a_rs[i] = el.mSymbol[i];
      v.push_back(rsexh::Encode(a_rs, code.mGf));
      // CRC
      int crc = 0;
      for (const auto& el : v.back()) {
         crc ^= el;
      }
      v.back().push_back(crc);
      //
   }
   // Channel
   const double ber = 0.001;
   for (int pos1 = 0; auto& el : v) {
      for (int pos2 = 0; auto& el2 : el) {
         for (int i=0; i<4; ++i) {
            const bool error = roll_error(ber);
            if (error)
               std::cout << "Make error: pos1: " << pos1 << ", pos2: " << pos2 << '\n';
            el2 ^= (static_cast<int>(error) << i);
         }
         pos2++;
      }
      pos1++;
   }
   // Decode
   std::cout << "Decode RS\n";
   for (int i = 0; auto& el : v) {
      // std::cout << "i: " << i << '\n';
      // CRC
      int crc = 0;
      for (const auto& el2 : el) {
         crc ^= el2;
      }
      if (crc)
         std::cout << "First crc: " << crc << '\n';
      auto c = rsexh::CalculateSyndrome(el, code.R, code.mGf);
      bool is_ok = true;
      for (const auto& el_c: c) {
         is_ok &= el_c == 0;
      }
      if (!is_ok) {
         std::cout << "First cyndrome check: " << "Failure" << '\n';
         if (auto it = code.mLut_1_errors.find(c); it != code.mLut_1_errors.end()) {
            std::cout << "Correction\n";
            const auto [pos, corrector] = it.operator*().second;
            const int channel_value = el.at(pos);
            el[pos] = code.mGf.Sub(channel_value - 1, corrector) + 1;
            // CRC
            crc = 0;
            for (const auto& el2 : el) {
               crc ^= el2;
            }
            if (crc)
               std::cout << "Second crc: " << crc << '\n';
            auto c = rsexh::CalculateSyndrome(el, code.R, code.mGf);
            is_ok = true;
            for (const auto& el_c: c) {
               is_ok &= el_c == 0;
            }
            std::cout << "Second cyndrome check: " << (is_ok ? "Ok" : "Failure") << '\n';
         }
      }
      if (is_ok) {
         auto a_dec = rsexh::Decode(el, code.R, code.mGf);
         for (int j=0; j<a_dec.size(); ++j) {
            a_received[i].mSymbol[j] = a_dec.at(j);
            a_received[i].mStatus = hamming::SymbolStatus::Normal;
         }
      } else {
         std::cout << "Erasured\n";
         for (int j=0; j<code.M2; ++j) {
            a_received[i].mSymbol[j] = -1;
            a_received[i].mStatus = hamming::SymbolStatus::Erased;
         }
      }
      i++;
   }
   std::cout << "Decode Hamming\n";
   int erased;
   const bool is_ok_hamming = code.mHammingCode.Decode(a_received, erased);
   std::cout << (is_ok_hamming ? "Ok" : "Failure") << ", were erased: " << erased << std::endl;
   std::cout << "Check input equality\n";
   bool is_equal = true;
   for (int i = 0; i < code.mHammingCode.K; ++i) {
      is_equal &= a.at(i) == a_received.at(i);
   }
   std::cout << (is_equal ? "Ok" : "Failure") << std::endl;
   // hamming::show_codeword(a_received, code, "Decoded input a_received: ");
}


int main( int argc, char* argv[] )
{
   measure_ber();
   // test_rs();
   return 0;
}
