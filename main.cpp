#include <iostream>
#include <random>
#include <cassert>
#include <set>
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

void test_ex_hamming_code(bool is_systematic) {
   std::cout << "Test Extended Hamming (default) code: " << (is_systematic ? "systematic" : "nonsystematic") << std::endl;;
   static constexpr int R2 = 6;  // Количество проверочных символов внешнего кода.
   static constexpr int M2 = 9;  // Количество внутренних символов внешнего кода.
   static hamming::HammingExtended< int, R2, M2 > mHammingCode;
   std::vector<std::set<int>> test_erasures = {{2, 5, 20}, {3, 7, 17}, {2, 3, 14}, {11, 14}, {1, 2, 9, 12}};
   mHammingCode.SwitchToSystematic(is_systematic);
   hamming::CodeWord<int, M2> a(mHammingCode.K);
   for (int erasure_round = 0; const auto& erasures : test_erasures) {
      std::cout << " ... round: " << erasure_round << "... ";
      erasure_round++;
      // Source
      for (auto& el : a) {
         el.mStatus = hamming::SymbolStatus::Normal;
         for (auto& symbol : el.mSymbol)
            symbol = roll_uint() & 15; // Полубайты.
            // symbol = 0; // Полубайты.
      }
      // Hamming encode
      // std::cout << "Hamming encode\n";
      auto s_h = mHammingCode.Encode(a);
      // hamming::show_codeword(s_h, mHammingCode.K, "Channel input v: ");
      for (int i = 0; auto& el : s_h) {
         if (erasures.contains(i)) {
            el.mStatus = hamming::SymbolStatus::Erased;
            for (auto& s : el.mSymbol)
               s = -1;
         }
         ++i;
      }
      // hamming::show_codeword(s_h, mHammingCode.K, "Channel output v: ");
      int erased;
      int was_changed_strategy;
      const bool is_ok_hamming = mHammingCode.Decode(s_h, erased, was_changed_strategy);
      bool is_equal = true;
      for (int i = 0; i < mHammingCode.K; ++i) {
         is_equal &= a.at(i) == s_h.at(i);
      }
      if (!is_equal) {
         std::cout << "Erased: " << erased << ", is ok hamming: " <<
                  is_ok_hamming << std::endl;
         std::cout << "Is equal: " << is_equal << std::endl;
         hamming::show_codeword(a, mHammingCode.K, "Input a: ");
         hamming::show_codeword(s_h, mHammingCode.K, "Decoded a: ");
      }
      std::cout << (is_equal ? "Ok." : "Failure.") << std::endl;
      if (!is_equal)
         break;
   }
}

void test_golay_code(bool is_systematic) {
   std::cout << "Test Golay code: " << (is_systematic ? "systematic" : "nonsystematic") << std::endl;
   static constexpr int R2 = 11; // Количество проверочных символов внешнего кода.
   static constexpr int M2 = 9;  // Количество внутренних символов внешнего кода.
   static hamming::HammingExtended< int, R2, M2 > mHammingCode{
      rsexh::Matrix<int>{ // Код Голея, циклический. Кодовое расстояние - 7.
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1}, 
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0}, 
            {0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0}, 
            {0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0}, 
            {0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0}, 
            {0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0}, 
            {0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0}, 
            {0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0}, 
            {0, 0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0}, 
            {0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
            {1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}
      }, 7
   };
   std::vector<std::set<int>> test_erasures = {{2, 5, 20}, {1, 6, 9, 12}, {3, 7, 17}, {2, 3, 14}, {0, 4, 13, 15, 16}, 
         {10, 11, 16, 17}, {4, 9, 10, 11, 14}, {0, 1, 6, 9, 11}, {0, 2, 5, 6, 8, 10}, {1, 3, 7, 19}, {0, 8, 9, 16, 21}};
   mHammingCode.SwitchToSystematic(is_systematic);
   hamming::CodeWord<int, M2> a(mHammingCode.K);
   for (int erasure_round = 0; const auto& erasures : test_erasures) {
      std::cout << " ... round: " << erasure_round << "... ";
      erasure_round++;
      // Source
      for (auto& el : a) {
         el.mStatus = hamming::SymbolStatus::Normal;
         for (auto& symbol : el.mSymbol)
            symbol = roll_uint() & 15; // Полубайты.
            // symbol = 0; // Полубайты.
      }
      // Encode
      auto s_h = mHammingCode.Encode(a);
      // hamming::show_codeword(s_h, mHammingCode.K, "Channel input v: ");
      for (int i = 0; auto& el : s_h) {
         if (erasures.contains(i)) {
            el.mStatus = hamming::SymbolStatus::Erased;
            for (auto& s : el.mSymbol)
               s = -1;
         }
         ++i;
      }
      // hamming::show_codeword(s_h, mHammingCode.K, "Channel output v: ");
      int erased;
      int was_changed_strategy;
      const bool is_ok_hamming = mHammingCode.Decode(s_h, erased, was_changed_strategy);
      bool is_equal = true;
      for (int i = 0; i < mHammingCode.K; ++i) {
         is_equal &= a.at(i) == s_h.at(i);
      }
      if (!is_equal) {
         std::cout << "Erased: " << erased << ", is ok hamming: " <<
                  is_ok_hamming << std::endl;
         std::cout << "Is equal: " << is_equal << std::endl;
         hamming::show_codeword(a, mHammingCode.K, "Input a: ");
         hamming::show_codeword(s_h, mHammingCode.K, "Decoded a: ");
      }
      std::cout << (is_equal ? "Ok." : "Failure.") << std::endl;
      if (!is_equal)
         break;
   }
}

void test_rs(int input) {
   rsexh::RsExh code;

   std::cout << "N: " << code.N << '\n';

   std::vector<int> a(code.K, input);
   auto s = rsexh::Encode(a, code.mGf);
   auto c = rsexh::CalculateSyndrome(s, code.R, code.mGf);

   auto a_dec = rsexh::Decode(s, code.R, code.mGf);

   rsexh::show_vector(a, "RS input:");
   rsexh::show_vector(s, "RS output:");
   rsexh::show_vector(c, "cyndrome:");
   rsexh::show_vector(a_dec, "RS decoded:");
}

void test_rs_correct_1(int position, int val) {
   static rsexh::RsExh code;
   position %= code.N;
   position = (code.N + position) % code.N;
   val = std::abs(val);
   val %= 16; // Полубайты; элементы поля Галуа GF(2^4).
   if (val == 0) {
      val += 1;
   }
   // std::cout << "\nRS 1-error correction test, error position: " << position << '\n';
   std::vector<int> a(code.K, 0);
   auto s = rsexh::Encode(a, code.mGf);
   // rsexh::show_vector(a, "RS input:");
   // rsexh::show_vector(s, "RS output:");
   for (int pos = 0; auto& el : s) {
      const bool error = (pos == position);
      if (error)
         el ^= val; // Полубайт.
      pos++;
   }
   // rsexh::show_vector(s, "Channel output:");
   auto c = rsexh::CalculateSyndrome(s, code.R, code.mGf);
   bool is_ok = true;
   for (const auto& el_c: c) {
      is_ok &= el_c == 0;
   }
   // rsexh::show_vector(c, "cyndrome 1:");
   if (!is_ok) {
      if (auto it = code.mLut_1_errors.find(c); it != code.mLut_1_errors.end()) {
         // std::cout << "Correction 1-error\n";
         const auto [pos, corrector_idx] = it.operator*().second;
         const int channel_value = s.at(pos);
         s[pos] = code.mGf.Sub(channel_value - 1, corrector_idx) + 1; // idx = value - 1 => value = idx + 1.
         c = rsexh::CalculateSyndrome(s, code.R, code.mGf);
         is_ok = true;
         for (const auto& el_c: c) {
            is_ok &= el_c == 0;
         }
         assert(is_ok);
      }
   }
   if (is_ok) {
      auto a_dec = rsexh::Decode(s, code.R, code.mGf);
      // rsexh::show_vector(s, "RS output:");
      // rsexh::show_vector(c, "cyndrome 2:");
      // rsexh::show_vector(a_dec, "RS decoded:");
      // std::cout << "Ok" << std::endl;
      assert(a_dec == a);
   } else {
      std::cout << "Failure 1-error correction." << std::endl;
      assert( 1 == 0 );
   }
}

void test_rs_correct_2(int position_1, int position_2, int val_1, int val_2) {
   static rsexh::RsExh code;
   val_1 = std::abs(val_1);
   val_2 = std::abs(val_2);
   position_1 %= code.N;
   position_2 %= code.N;
   position_1 = (code.N + position_1) % code.N;
   position_2 = (code.N + position_2) % code.N;
   if (position_1 == position_2)
      return;
   val_1 %= 16; // Полубайты; элементы поля Галуа GF(2^4).
   val_2 %= 16;
   if (val_1 == 0) {
      val_1 += 1;
   }
   if (val_2 == 0) {
      val_2 += 1;
   }
   // std::cout << "\nRS 2-error correction test, positions: " << position_1 << ", " << position_2 << ", values: " << val_1 << ", " << val_2 << '\n';
   std::vector<int> a(code.K, 0);
   auto s = rsexh::Encode(a, code.mGf);
   // rsexh::show_vector(a, "RS input:");
   // rsexh::show_vector(s, "RS output:");
   for (int pos = 0; auto& el : s) {
      const bool error1 = (pos == position_1);
      const bool error2 = (pos == position_2);
      if (error1)
         el ^= val_1;
      if (error2)
         el ^= val_2;
      pos++;
   }
   // rsexh::show_vector(s, "Channel output:");
   auto c = rsexh::CalculateSyndrome(s, code.R, code.mGf);
   bool is_ok = true;
   for (const auto& el_c: c) {
      is_ok &= el_c == 0;
   }
   // rsexh::show_vector(c, "cyndrome 1:");
   if (!is_ok) {
      for (int k = 0; k < code.N - 1; k++) {
         if (auto it = code.mLut_2_errors.find(c); it != code.mLut_2_errors.end()) {
            const auto [pos_2nd, corrector_indices] = it.operator*().second;
            const int idx_1 = k;
            const int idx_2 = pos_2nd + k;
            // std::cout << "Correction 2-error: " << pos_2nd << ", " << k << std::endl;
            const int channel_value_1 = s.at(idx_1);
            const int channel_value_2 = s.at(idx_2);
            const auto [corrector_idx_1, corrector_idx_2] = corrector_indices;
            // std::cout << "Correction indices: " << corrector_idx_1 << ", " << corrector_idx_2 << std::endl;
            s[idx_1] = code.mGf.Sub(channel_value_1 - 1, corrector_idx_1) + 1;
            s[idx_2] = code.mGf.Sub(channel_value_2 - 1, corrector_idx_2) + 1; 
            c = rsexh::CalculateSyndrome(s, code.R, code.mGf);
            is_ok = true;
            for (const auto& el_c: c) {
               is_ok &= el_c == 0;
            }
            assert(is_ok);
            break;
         }
         rsexh::ShiftLeftSyndrome<code.p, code.q>(c); // Сдвиг - имеется ввиду сдвиг соответствующего вектора ошибки.
      }
   }
   if (is_ok) {
      auto a_dec = rsexh::Decode(s, code.R, code.mGf);
      // rsexh::show_vector(s, "Corrected channel output:");
      // rsexh::show_vector(c, "cyndrome 2:");
      // rsexh::show_vector(a_dec, "RS decoded:");
      assert(a_dec == a);
      // std::cout << "Ok" << std::endl;
   } else {
      std::cout << "Failure 2-error correction." << std::endl;
      assert( 1 == 0 );
   }
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
      // hamming::show_codeword(s_h, code.mHammingCode.K, "Hamming output s: ");
      // RS encode
      // std::cout << "RS encode\n";
      v.clear();
      for (const auto& el : s_h) {
         assert(code.K == el.mSymbol.size());
         for (int i = 0; i < code.K; ++i)
            a_rs[i] = el.mSymbol[i];
         v.push_back(rsexh::Encode(a_rs, code.mGf));
      }
      // rsexh::show_matrix(v, "RS outputs: ");
      // Channel
      std::vector<int> error_q; // Кратности ошибки.
      for (int pos1 = 0; auto& el : v) {
         error_q.push_back(0);
         for (int pos2 = 0; auto& el2 : el) {
            bool was_error = false;
            for (int i=0; i<4; ++i) {
               const bool error = roll_error(ber);
               was_error |= error;
               el2 ^= (static_cast<int>(error) << i); // Полубайт.
            }
            error_q.back() += was_error;
            pos2++;
         }
         pos1++;
      }
      // Decode
      std::vector<int> was_1_error_correction(v.size());
      std::vector<int> was_2_error_correction(v.size());
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
               was_1_error_correction[i] = 1;
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
                  was_2_error_correction[i] = 1;
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
      // hamming::show_codeword(a_received, code.mHammingCode.K, "Received Hamming input v: ");
      [[maybe_unused]] int was_changed_strategy;
      const bool is_ok_hamming = code.mHammingCode.Decode(a_received, erased, was_changed_strategy);
      // if (erased) {
      //    std::cout << (is_ok_hamming ? "Ok" : "Failure") << ", were erased: " << erased << std::endl;
      // }
      // std::cout << "Check input equality\n";
      bool is_equal = true;
      for (int i = 0; i < code.mHammingCode.K; ++i) {
         is_equal &= a.at(i) == a_received.at(i);
      }
      // if (was_changed_strategy) {
      //    std::cout << "is equal after strategy changing: " << (is_equal ? "yes" : "no") << std::endl;
      // }
      if (!is_equal) {
         const int code_distance = code.mHammingCode.D;
         int bad_correction = 0; // Была неисправимая ошибка.
         for (int pos = 0; const auto el : was_1_error_correction) {
            if (el && error_q.at(pos) != 1) {
               bad_correction = 1;
               break;
            }
            pos++;
         }
         for (int pos = 0; const auto el : was_2_error_correction) {
            if (el && error_q.at(pos) != 2) {
               bad_correction = 1;
               break;
            }
            pos++;
         }
         if ((erased < code_distance) && !bad_correction) { // Если ошибка исправима, то стирание гарантированной кратности должно быть восстановлено.
            std::cout << "Failure: erased: " << erased << ", is ok hamming: " << is_ok_hamming << std::endl;
            rsexh::show_vector(was_1_error_correction, "1-error corrections");
            rsexh::show_vector(was_2_error_correction, "2-error corrections");
            rsexh::show_vector(error_q, "Channel errors (q)");
            hamming::show_codeword(a, code.mHammingCode.K, "Input a: ");
            hamming::show_codeword(a_received, code.mHammingCode.K, "Decoded a: ");
            rsexh::show_matrix(code.mHammingCode.mErasureSubmatrix, "Selected matrix: ");
            return -1.;
         }
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
   test_golay_code(true);
   test_golay_code(false);

   test_ex_hamming_code(true);
   test_ex_hamming_code(false);

   // Channel BER : Decoder BER
   
   // Case A.
   // RS (15, 10, 6) in mode 1- and 2-error correction.
   // Default Extended Hamming code (32, 26) with the code distance 4. Total R = 46%.
   // 0.002 : 
   // 0.005 : 
   // 0.010 : 
   // 0.015 : 
   // 0.020 : 
   // 0.025 : 
   // 0.030 : 
   // 0.100 : 
   // 0.200 : 
   
   // Case B.
   // RS (15, 9, 7) in mode 1- and 2-error correction.
   // Set the Golay code (23, 12) with the code distance 7. Total R = 69%.
   // 0.005 : 
   // 0.010 : 
   // 0.015 : 6.0e-7
   // 0.020 : 1.5e-5
   // 0.025 : 2.8e-4
   // 0.030 : 0.0041
   // 0.100 : 0.465

   const double ber = 0.02;
   double output_ber = 0;
   for (double counter = 1;; counter++) {
      double prev_ber = output_ber;
      auto out_ber = measure_ber(ber, 10000);
      if (out_ber == -1.) {
         return -1;
      }
      output_ber += (out_ber - output_ber) / counter;
      const double rel_error = output_ber != 0. ? std::abs(prev_ber - output_ber) / output_ber : 1.;
      std::cout << "decoder BER: " << output_ber << "\tcounter: " << counter << "\tchannel BER: " << ber << "\t(sample ber = " << out_ber << ")" << std::endl;
      if (out_ber > 0 && rel_error < 1.e-4) {
         break;
      }
   }
   std::cout << "Decoder BER: " << output_ber << std::endl;

   // test_rs(0);
   // test_rs(15);

   // for (;;) {
   //    test_rs_correct_1(roll_uint(), roll_uint());
   //    test_rs_correct_2(roll_uint(), roll_uint(), roll_uint(), roll_uint());
   // }
   return 0;
}
