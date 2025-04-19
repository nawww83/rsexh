#pragma once

#include <cstdint>
#include <cassert>
#include <vector>

namespace lfsr8
{

class LFSR
{
   using STATE = std::vector< int >;
   using SAMPLE = int;

public:
   explicit LFSR( const STATE& K, int p = 2 )
      : m_K( K )
   {
      mP = p;
      mQ = K.size();
      m_state.resize( K.size() );
      m_calculate_inverse_of_K();
   };

   void set_state( const STATE& st )
   {
      m_state = st;
   }

   void set_unit_state()
   {
      const auto old_size = m_state.size();
      m_state.clear();
      m_state.resize( old_size, 0 );
      m_state[ 0 ] = 1;
   }

   void set_K( STATE K )
   {
      m_K = K;
      mQ = K.size();
      m_calculate_inverse_of_K();
   }

   /**
    * @brief Сделать шаг вперед (один такт генератора).
    * @param input Входной символ (по модулю p), который подается
    * на вход генератора.
    */
   void next( SAMPLE input = 0 )
   {

      const SAMPLE m_v = m_state[ mQ - 1 ];
      for( int i = mQ - 1; i > 0; i-- )
      {
         m_state[ i ] = ( m_state[ i - 1 ] + m_v * m_K[ i ] ) % SAMPLE( mP );
      }
      m_state[ 0 ] = ( input + m_v * m_K[ 0 ] ) % SAMPLE( mP );
   }

   /**
    * @brief Сделать шаг назад (один такт генератора). Обратно к next(input).
    * @param input Входной символ (по модулю p), который подается
    * на вход генератора.
    */
   void back( SAMPLE input = 0 )
   {
      const SAMPLE m_v = ( m_inv_K0 * ( m_state[ 0 ] - input + SAMPLE( mP ) ) ) % SAMPLE( mP );
      for( int i = 0; i < mQ - 1; i++ )
      {
         m_state[ i ] = ( m_state[ i + 1 ] - m_v * m_K[ i + 1 ] + SAMPLE( mP ) * SAMPLE( mP ) ) % SAMPLE( mP );
      }
      m_state[ mQ - 1 ] = m_v;
   }

   /**
    * @brief Является ли заданное состояние текущим состоянием генератора.
    * @param st Заданное состояние.
    * @return Да/нет.
    */
   bool is_state( const STATE& st ) const
   {
      return st == m_state;
   }

   auto get_state() const
   {
      return m_state;
   }

   auto get_cell( int idx ) const
   {
      return m_state[ idx ];
   }

private:
   STATE m_state{};
   STATE m_K{};
   SAMPLE m_inv_K0{};
   int mP;
   int mQ;

   /**
    * @brief Вычисляется обратный (по умножению) коэффициент.
    */
   void m_calculate_inverse_of_K()
   {
      const auto x = m_K[ 0 ];
      assert( x != 0 );
      assert( mP > 1 );
      m_inv_K0 = 1;
      for( ;; )
      {
         const auto modulo_p = x * m_inv_K0 % static_cast< SAMPLE >( mP );
         if( modulo_p == SAMPLE( 1 ) )
         {
            break;
         }
         m_inv_K0++;
      }
   }
};

} // namespace lfsr8
