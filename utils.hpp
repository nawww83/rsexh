#pragma once


namespace utils
{

/**
 * Степень 2.
 */
template <typename T>
inline constexpr T power2( int x )
{
   if ( x >= 0 )
      return T( 1 ) << x;
   else
      return T(1);
}

template <>
inline constexpr int power2( int x )
{
   if ( x >= 0 )
      return 1 << x;
   else
      return 1;
}

/**
 * Степень p.
 */
template <typename T>
inline constexpr T power(T p, int x )
{
    T result = 1;
    while ( x > 0 ) {
        result *= p;
        x--;
    }
    return result;
}

} // namespace utils