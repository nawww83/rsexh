#pragma once

#include <unordered_map>
#include <map>
#include <vector>

namespace gf {

/**
 * @brief Состояния в виде вектора элементов поля GF(p).
 * Это могут быть как коэффициенты примитивного полинома, порождающего поле
 * Галуа GF(p^q), так и элементы самого поля GF(p^q) в векторной форме.
 */
struct State
{
   explicit State() = default;
   State( const State& ) = default;
   State( State&& ) = default;
   State& operator=( const State& other ) = default;
   State& operator=( State&& other ) = default;

   State( int p, int n )
      : mP{ p }
   {
      mState.clear();
      mState.resize( n, 0 );
   }
   State( int p, const std::vector< int >& state )
      : mP{ p }
      , mState{ state }
   {
   }

   int mP = 0;
   std::vector< int > mState;

   State operator+( const State& other ) const;

   State operator-( const State& other ) const;

   bool operator==( const State& other ) const;

   auto Size() const
   {
      return mState.size();
   }

   void Resize( std::size_t size )
   {
      mState.resize( size );
   }
};

struct KeyHasher1
{
   std::size_t operator()( const State& st ) const
   {
      using std::hash;
      auto tmp = hash< int >()( st.mP );
      for( auto& el : st.mState )
      {
         tmp ^= hash< int >()( el );
      }
      return tmp;
   }
};

struct KeyHasher2
{
   std::size_t operator()( const std::vector< int >& st ) const
   {
      using std::hash;
      auto tmp = 0;
      for( auto& el : st )
      {
         tmp ^= hash< int >()( el );
      }
      return tmp;
   }
};

/**
 * @brief Таблица соответствия между множествами состояний и индексов.
 * Некоторый индекс - это степень элемента "альфа" поля Галуа GF(p^q).
 */
class GFLUT
{
public:
   explicit GFLUT() = default;

   /**
    * Конструктор, принимающий на вход порождающий полином.
    */
   GFLUT( const State& g_poly );

   int Index( const State& st ) const;

   State Element( int idx ) const;

   auto Size() const;

   bool PolyIsGood() const;

   /**
    * Формирует и возвращает упорядоченную по индексу таблицу 
    * соответствий между индексами и состояниями.
    */
   std::map< int, State > OrderedLut() const;

private:
   /**
    * Прямая таблица соответствий.
    */
   std::unordered_map< State, int, KeyHasher1 > mLUT;

   /**
    * Обратная таблица соответствий.
    */
   std::unordered_map< int, State > mInvLUT;
   
   /**
    * Флаг, определяющий примитивен ли порождающий полином.
    */
   bool mPolyIsGood = false;

    /**
    * Заполняет по порождающему полиному таблицы соответствий.
    */
   std::size_t FillLut( const State& g_poly );
};

/**
 * Класс для арифметических манипуляций с полем Галуа GF(p^q).
 */
template< int p, int q >
class GF
{
public:
   explicit GF() = default;
   
   /**
    * Конструктор с переданной ссылкой на уже сформированную таблицу соответствия.
    */
   GF( std::reference_wrapper< GFLUT > lut );
   
   State Add( const State& lhs, const State& rhs ) const;
   int Add( const int idx1, const int idx2 ) const;
   
   State Sub( const State& lhs, const State& rhs ) const;
   int Sub( const int idx1, const int idx2 ) const;
   
   State Mult( const State& lhs, const State& rhs ) const;
   int Mult( const int idx1, const int idx2 ) const;
   
   int GetIndex( const State& st ) const;
   State GetElement( int idx ) const;

private:   
   std::reference_wrapper< GFLUT > mLUT;
};

}