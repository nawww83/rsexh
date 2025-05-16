#include <iostream>
#include <random>
#include <cassert>
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

void test_rs_correct_1(int position, int val, int pos_2 = -1, int val_2 = -1) {
   rsexh::RsExh code;
   position = (code.N + position) % code.N;
   std::cout << "RS 1-error correction test, error position: " << position << '\n';
   std::vector<int> a(code.K, 14);
   auto s = rsexh::Encode(a, code.mGf);
   rsexh::show_vector(a, "RS input:");
   rsexh::show_vector(s, "RS output:");
   for (int pos = 0; auto& el : s) {
      const bool error = (pos == position);
      const bool error_2 = (pos == pos_2);
      if (error)
         el ^= val; // Полубайт.
      if (error_2)
         el ^= val_2;
      pos++;
   }
   auto c = rsexh::CalculateSyndrome(s, code.R, code.mGf);
   rsexh::show_vector(c, "cyndrome:");
   if (auto it = code.mLut_1_errors.find(c); it != code.mLut_1_errors.end()) {
      std::cout << "Correction 1-error\n";
      const auto [pos, corrector_idx] = it.operator*().second;
      const int channel_value = s.at(pos);
      s[pos] = code.mGf.Sub(channel_value - 1, corrector_idx) + 1; // idx = value - 1 => value = idx + 1.
      c = rsexh::CalculateSyndrome(s, code.R, code.mGf);
      bool is_ok = true;
      for (const auto& el_c: c) {
         is_ok &= el_c == 0;
      }
      assert(is_ok);
   }
   auto a_dec = rsexh::Decode(s, code.R, code.mGf);
   rsexh::show_vector(s, "RS output:");
   rsexh::show_vector(c, "cyndrome:");
   rsexh::show_vector(a_dec, "RS decoded:");
   std::cout << (a_dec == a ? "Corrected" : "Failure") << std::endl;
   assert(a_dec == a);
}

void test_rs_correct_2(int position_1, int position_2, int val_1, int val_2, int position_3 = -1, int val_3 = -1) {
   rsexh::RsExh code;
   position_1 = (code.N + position_1) % code.N;
   position_2 = (code.N + position_2) % code.N;
   std::cout << "RS 2-error correction test, error positions: " << position_1 << ", " << position_2 << '\n';
   std::vector<int> a(code.K, 14);
   auto s = rsexh::Encode(a, code.mGf);
   rsexh::show_vector(a, "RS input:");
   rsexh::show_vector(s, "RS output:");
   for (int pos = 0; auto& el : s) {
      const bool error1 = (pos == position_1);
      const bool error2 = (pos == position_2);
      const bool error3 = (pos == position_3);
      if (error1)
         el ^= val_1;
      if (error2)
         el ^= val_2;
      if (error3)
         el ^= val_3;
      pos++;
   }
   rsexh::show_vector(s, "Channel output:");
   auto c = rsexh::CalculateSyndrome(s, code.R, code.mGf);
   rsexh::show_vector(c, "cyndrome:");
   for (int k = 0; k < code.N - 1; k++) {
      if (auto it = code.mLut_2_errors.find(c); it != code.mLut_2_errors.end()) {
         const auto [pos_2nd, corrector_indices] = it.operator*().second;
         const int idx_1 = k;
         const int idx_2 = pos_2nd + k;
         std::cout << "Correction 2-error: " << pos_2nd << ", " << k << std::endl;
         const int channel_value_1 = s.at(idx_1);
         const int channel_value_2 = s.at(idx_2);
         const auto [corrector_idx_1, corrector_idx_2] = corrector_indices;
         std::cout << "Correction indices: " << corrector_idx_1 << ", " << corrector_idx_2 << std::endl;
         s[idx_1] = code.mGf.Sub(channel_value_1 - 1, corrector_idx_1) + 1;
         s[idx_2] = code.mGf.Sub(channel_value_2 - 1, corrector_idx_2) + 1; 
         c = rsexh::CalculateSyndrome(s, code.R, code.mGf);
         bool is_ok = true;
         for (const auto& el_c: c) {
            is_ok &= el_c == 0;
         }
         assert(is_ok);
         break;
      }
      rsexh::ShiftLeftSyndrome<code.p, code.q>(c); // Сдвиг - имеется ввиду сдвиг соответствующего вектора ошибки.
   }
   auto a_dec = rsexh::Decode(s, code.R, code.mGf);
   rsexh::show_vector(s, "Corrected channel output:");
   rsexh::show_vector(c, "cyndrome:");
   rsexh::show_vector(a_dec, "RS decoded:");
   std::cout << (a_dec == a ? "Corrected" : "Failure") << std::endl;
   assert(a_dec == a);
}

double measure_ber(double ber, int factor) {
   static rsexh::RsExh code;
   code.mHammingCode.SwitchToSystematic(false);
   // std::cout << "N: " << code.N << '\n';
   hamming::CodeWord<int, code.M2> a(code.mHammingCode.K);
   hamming::CodeWord<int, code.M2> a_received(code.mHammingCode.N);
   long long bits_transmitted = 0;
   long long bits_corrupted = 0;
   double output_ber = 0.;
   for (int f = 0; f < factor; f++) {
      a_received.resize(code.mHammingCode.N);
      // Source
      for (auto& el : a) {
         el.mStatus = hamming::SymbolStatus::Normal;
         for (auto& symbol : el.mSymbol)
            symbol = roll_uint() & 15; // Полубайты.
      }
      bits_transmitted += code.mHammingCode.K * code.M2 * 4;
      // hamming::show_codeword(a, code, "Input a: ");
      // Hamming encode
      // std::cout << "Hamming encode\n";
      auto s_h = code.mHammingCode.Encode(a);
      std::vector<int> a_rs(code.K);
      std::vector<std::vector<int>> v;
      // RS encode
      // std::cout << "RS encode\n";
      v.clear();
      for (const auto& el : s_h) {
         assert(code.K == el.mSymbol.size());
         for (int i = 0; i < code.K; ++i)
            a_rs[i] = el.mSymbol[i];
         v.push_back(rsexh::Encode(a_rs, code.mGf));
      }
      // Channel
      // std::vector<int> error_q; // Кратности ошибки.
      for (int pos1 = 0; auto& el : v) {
         // error_q.push_back(0);
         for (int pos2 = 0; auto& el2 : el) {
            bool was_error = false;
            for (int i=0; i<4; ++i) {
               const bool error = roll_error(ber);
               // was_error |= error;
               el2 ^= (static_cast<int>(error) << i); // Полубайт.
            }
            // error_q.back() += was_error;
            pos2++;
         }
         pos1++;
      }
      // Decode
      for (int i = 0; auto& el : v) {
         a_received[i].mStatus = hamming::SymbolStatus::Uninitialized;
         // std::cout << "Calculate cyndrome 1\n";
         auto c = rsexh::CalculateSyndrome(el, code.R, code.mGf);
         bool is_ok = true;
         for (const auto& el_c: c) {
            is_ok &= el_c == 0;
         }
         if (!is_ok) {
            // std::cout << "First cyndrome check: " << "Failure" << '\n';
            if (auto it = code.mLut_1_errors.find(c); it != code.mLut_1_errors.end()) {
               const auto [pos, corrector_idx] = it.operator*().second;
               const int channel_value = el.at(pos);
               el[pos] = code.mGf.Sub(channel_value - 1, corrector_idx) + 1; // idx = value - 1 => value = idx + 1.
               is_ok = true;
            }
         }
         if (is_ok) {
            auto a_dec = rsexh::Decode(el, code.R, code.mGf);
            for (int j=0; j<a_dec.size(); ++j) {
               a_received[i].mSymbol[j] = a_dec.at(j);
            }
            a_received[i].mStatus = hamming::SymbolStatus::Normal;
         } else {
            for (int j=0; j<code.M2; ++j) {
               a_received[i].mSymbol[j] = -1;
            }
            a_received[i].mStatus = hamming::SymbolStatus::Erased;
         }
         i++;
      }
      for (int i = 0; auto& el : v) {
         if (a_received.at(i).mStatus == hamming::SymbolStatus::Normal) {
            i++;
            continue;
         }
         // std::cout << "Calculate cyndrome 2\n";
         auto c = rsexh::CalculateSyndrome(el, code.R, code.mGf);
         bool is_ok = true;
         for (const auto& el_c: c) {
            is_ok &= el_c == 0;
         }
         if (!is_ok) {
            // std::cout << "First cyndrome check: " << "Failure" << '\n';
            // 2-кратные ошибки.
            for (int k = 0; k < code.N - 1; k++) {
               if (auto it = code.mLut_2_errors.find(c); it != code.mLut_2_errors.end()) {
                  const auto [pos_2nd, corrector_indices] = it.operator*().second;
                  const int idx_1 = k;
                  const int idx_2 = pos_2nd + k;
                  if (idx_2 < el.size()) {
                     const int channel_value_1 = el.at(idx_1);
                     const int channel_value_2 = el.at(idx_2);
                     const auto [corrector_idx_1, corrector_idx_2] = corrector_indices;
                     el[idx_1] = code.mGf.Sub(channel_value_1 - 1, corrector_idx_1) + 1;
                     el[idx_2] = code.mGf.Sub(channel_value_2 - 1, corrector_idx_2) + 1;
                  }
                  is_ok = true;
                  break;
               }
               rsexh::ShiftLeftSyndrome<code.p, code.q>(c); // Сдвиг - имеется ввиду сдвиг соответствующего вектора ошибки.
            }
         }
         if (is_ok) {
            auto a_dec = rsexh::Decode(el, code.R, code.mGf);
            for (int j=0; j<a_dec.size(); ++j) {
               a_received[i].mSymbol[j] = a_dec.at(j);
            }
            a_received[i].mStatus = hamming::SymbolStatus::Normal;
         } else {
            for (int j=0; j<code.M2; ++j) {
               a_received[i].mSymbol[j] = -1;
            }
            a_received[i].mStatus = hamming::SymbolStatus::Erased;
         }
         i++;
      }
      // std::cout << "Decode Hamming: " << a_received.size() << std::endl;
      int erased;
      const bool is_ok_hamming = code.mHammingCode.Decode(a_received, erased);
      // if (erased) {
         // std::cout << (is_ok_hamming ? "Ok" : "Failure") << ", were erased: " << erased << std::endl;
      // }
      // std::cout << "Check input equality\n";
      bool is_equal = true;
      for (int i = 0; i < code.mHammingCode.K; ++i) {
         is_equal &= a.at(i) == a_received.at(i);
      }
      if (!is_equal) {
         for (int i = 0; i < code.mHammingCode.K; ++i) {
            if (a_received.at(i).mStatus != hamming::SymbolStatus::Normal) {
               bits_corrupted += ((code.M2*4) / 2); // Все полубайты повреждены: стерты, эквивалентная вероятность ошибки 0.5.
               continue;
            }
            for (int j = 0; j < code.M2; ++j) {
               int error = a.at(i).mSymbol.at(j) ^ a_received.at(i).mSymbol.at(j);
               bits_corrupted += (error & 1);
               error >>= 1;
               bits_corrupted += (error & 1);
               error >>= 1;
               bits_corrupted += (error & 1);
               error >>= 1;
               bits_corrupted += (error & 1);
            }
         }
      }
   }
   output_ber = (1. * bits_corrupted) / bits_transmitted;
   return output_ber;
   // std::cout << (is_equal ? "Ok" : "Failure") << std::endl;
   // hamming::show_codeword(a_received, code, "Decoded input a_received: ");
}


int main( int argc, char* argv[] )
{
   // Channel BER : decoder BER, RS (15,10,6) in mode 1-error and 2-error correction. Extended Hamming code (32, 26) with distance 4.
   // 0.002 : 1.0e-7
   // 0.005 : 2.3e-6
   // 0.010 : 0.00023
   // 0.015 : 0.0047
   // 0.020 : 0.026
   // 0.025 : 0.070
   // 0.030 : 0.117
   // 0.100 : 0.467
   // 0.200 : 0.5
   
   // Set the Golay code (23, 12) with code distance 7. Approx., BER 10 times less.
   // 0.010 : 1.7e-5
   // 0.015 : 0.00044

   const double ber = 0.015;
   double output_ber = 0;
   for (double counter = 1;; counter++) {
      double prev_ber = output_ber;
      auto out_ber = measure_ber(ber, 10000);
      if (out_ber == -1.) {
         return -1;
      }
      output_ber += (out_ber - output_ber) / counter;
      const double rel_error = output_ber != 0. ? std::abs(prev_ber - output_ber) / output_ber : 1.;
      std::cout << "decoder BER: " << output_ber << ", counter: " << counter << ", channel BER: " << ber << std::endl;
      if (out_ber > 0 && rel_error < 1.e-4) {
         break;
      }
   }
   std::cout << "Decoder BER: " << output_ber << std::endl;

   // test_rs();
   
   // test_rs_correct_1(0);
   // test_rs_correct_1(1);
   // test_rs_correct_1(-1);
   // 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 8
   // test_rs_correct_1(14, 2, 15, 8);

   // test_rs_correct_2(0, 1, 15, 15);
   // test_rs_correct_2(1, 2, 15, 15);
   // test_rs_correct_2(-1, -2, 15, 15);
   // 0, 0, 1, 0, 0, 0, 0, 0, 0, 4, 0, 1, 0, 0, 0
   // test_rs_correct_2(2, 9, 4, 1, 11, 1);
   return 0;
}
