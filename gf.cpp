#include <cmath>
#include <cassert>

#include "lfsr.hpp"
#include "gf.hpp"

namespace gf {

GFLUT::GFLUT( const State& g_poly )
{
   const auto Order = FillLut( g_poly );
   if( mInvLUT.size() < Order || mLUT.size() < Order )
   {
      mPolyIsGood = false;    
      return;
   }
   mPolyIsGood = true;
}

int GFLUT::Index( const State& st ) const
{
   return mLUT.at( st );
}

State GFLUT::Element( int idx ) const
{
   return mInvLUT.at( idx );
}

auto GFLUT::Size() const
{
   return std::min( mLUT.size(), mInvLUT.size() );
}

bool GFLUT::PolyIsGood() const
{
   return mPolyIsGood;
}

std::map< int, State > GFLUT::OrderedLut() const
{
   std::map< int, State > result;
   for( const auto& [ idx, state ] : mInvLUT )
   {
      result[ idx ] = state;
   }
   return result;
}

std::size_t GFLUT::FillLut( const State& g_poly )
{
   const int q = g_poly.Size();
   const int p = g_poly.mP;
   const std::size_t Order = std::pow( p, q );
   const auto& poly = g_poly.mState;
   lfsr8::LFSR gen{ poly, p };
   gen.set_unit_state();
   const State unit_state{ p, gen.get_state() };
   mLUT.clear();
   mInvLUT.clear();
   const State zero_st( p, q );
   mInvLUT[ -1 ] = zero_st;
   mLUT[ mInvLUT[ -1 ] ] = -1;
   mInvLUT[ 0 ] = { unit_state };
   mLUT[ mInvLUT[ 0 ] ] = 0;
   const auto& st_1{ unit_state.mState };
   for( int idx = 1;; idx++ )
   {
      gen.next();
      if( gen.is_state( st_1 ) )
      {
         break;
      }
      mInvLUT[ idx ] = { p, gen.get_state() };
      mLUT[ mInvLUT[ idx ] ] = idx;
   }
   return Order;
}

template< int p, int q >
GF< p, q >::GF( std::reference_wrapper< GFLUT > lut )
   : mLUT{ lut }
{
}

template< int p, int q >
int GF< p, q >::GetIndex( const State& st ) const
{
   return mLUT.get().Index( st );
}

template< int p, int q >
State GF< p, q >::GetElement( int idx ) const
{
   return mLUT.get().Element( idx );
}

template< int p, int q >
State GF< p, q >::Add( const State& lhs, const State& rhs ) const
{
   return GetElement( GetIndex( lhs + rhs ) );
}

template< int p, int q >
int GF< p, q >::Add( const int idx1, const int idx2 ) const
{
   const int N = mLUT.get().Size() - 1;
   return GetIndex( GetElement( idx1 % N ) + GetElement( idx2 % N ) );
}

template< int p, int q >
State GF< p, q >::Sub( const State& lhs, const State& rhs ) const
{
   return GetElement( GetIndex( lhs - rhs ) );
}

template< int p, int q >
int GF< p, q >::Sub( const int idx1, const int idx2 ) const
{
   const int N = mLUT.get().Size() - 1;
   return GetIndex( GetElement( idx1 % N ) - GetElement( idx2 % N ) );
}

template< int p, int q >
State GF< p, q >::Mult( const State& lhs, const State& rhs ) const
{
   const int N = mLUT.get().Size() - 1;
   const int idx1 = GetIndex( lhs );
   const int idx2 = GetIndex( rhs );
   return GetElement( idx1 >= 0 && idx2 >= 0 ? ( idx1 + idx2 ) % N : -1 );
}

template< int p, int q >
int GF< p, q >::Mult( const int idx1, const int idx2 ) const
{
   const int N = mLUT.get().Size() - 1;
   return idx1 >= 0 && idx2 >= 0 ? ( idx1 + idx2 ) % N : -1;
}

State State::operator+( const State& other ) const
{
   if( other.Size() != this->Size() )
   {
      return State{};
   }
   State result{ *this };
   for( std::size_t i = 0; i < result.Size(); ++i )
   {
      result.mState[ i ] += other.mState[ i ];
      result.mState[ i ] %= mP;
   }
   return result;
}

State State::operator-( const State& other ) const
{
   if( other.Size() != this->Size() )
   {
      return State{};
   }
   State result{ *this };
   for( std::size_t i = 0; i < result.Size(); ++i )
   {
      result.mState[ i ] -= other.mState[ i ];
      result.mState[ i ] += mP;
      result.mState[ i ] %= mP;
   }
   return result;
}

bool State::operator==( const State& other ) const
{
   return mP == other.mP && mState == other.mState;
}

template class GF< 2, 4 >;

}