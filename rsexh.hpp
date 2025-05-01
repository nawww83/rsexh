#pragma once

#include <cmath> // std::pow
#include <string> // std::string
#include <cassert> // assert
#include <iostream> // std::cout
#include "gf.hpp"
#include "hamming.hpp"
#include "utils.hpp" // power2

namespace rsexh {
    template <typename T>
    using Vector = std::vector<T>;

    template <typename T>
    using Matrix = Vector<Vector<T>>;

    /**
     * Показать вектор элементов.
     */
    template <typename T>
    static void show_vector(const Vector<T>& v, const std::string& title)
    {
        std::cout << title << '\n';    
        for (const auto& el : v) {
            std::cout << el << ", ";
        }
        std::cout << std::endl;    
    }

    /**
     * Показать матрицу элементов.
     */
    template <typename T>
    inline void show_matrix(const Matrix<T>& M, const std::string& title)
    {
        std::cout << title << '\n';
        for (const auto& row : M) {
            for (const auto& el : row) {
                std::cout << el << ", ";
            }
            std::cout << std::endl;
        }
    }

    /**
     * Сформировать проверочную матрицу кода Рида-Соломона (несистематический код).
     */
    template< int p, int q >
    inline Matrix<int> GetParityCheck( int R )
    {
        const int N = std::pow( p, q ) - 1;
        Matrix<int> result;
        for( int i = 0; i < R; ++i ) // Строки H.
        {
            int idx = 0;
            result.push_back( {} );
            for( int j = 0; j < N; ++j ) // Столбцы H.
            {
                result[ i ].push_back( idx );
                idx += ( i + 1 ); // Индекс элемена матрицы H.
                idx %= N;
            }
        }
        return result;
    }

    /**
     * Вычислить синдром по принятому вектору (несистематический код).
     */
    template< int p, int q >
    inline std::vector< int > CalculateSyndrome( const std::vector< int >& v, int R, const gf::GF< p, q >& gf )
    {   
        const int N = std::pow( p, q ) - 1;
        std::vector< int > result;
        for( int i = 0; i < R; ++i ) // Строки H.
        {
            int idx = 0;
            int result_idx = -1;
            for( int j = 0; j < N; ++j ) // Столбцы H.
            {
                const int mult_idx = gf.Mult( v.at( j ) - 1, idx );
                result_idx = gf.Add( mult_idx, result_idx );
                idx += ( i + 1 ); // Индекс элемента матрицы H.
                idx %= N;
            }
            result.push_back( result_idx + 1 );
        }
        return result;
    }

    /**
     * Несистематическое кодирование кодом Рида-Соломона.
     */
    template< int p, int q >
    inline std::vector< int > Encode( const std::vector< int >& a, const gf::GF< p, q >& gf)
    {
        const int N = std::pow( p, q ) - 1;
        std::vector< int > a_padded;
        for (auto el : a) {
            a_padded.push_back(el - 1);
        }
        std::vector< int > result;
        // Дополнение нулями справа. Минус единица - это индекс, которому соответствует нуль-элемент поля.
        const int K = a.size();
        const int R = N - K;
        for (int i = 0; i < R; ++i) {
            a_padded.push_back(-1);
        }
        // s = a'F, где F - квадратная матрица наподобие матрицы БПФ.
        // a' - дополненный нулями информационный вектор a.
        for (int i = 0; i < N; ++i) { // По столбцам матрицы F.
            const int step = i; // Шаг степени элементов "альфа".
            int idx = 0;
            int result_idx = -1;
            for (int j = 0; j < N; ++j) { // По строкам i-го столбца.
                const int mult_idx = gf.Mult( a_padded.at( j ), idx );
                result_idx = gf.Add( mult_idx, result_idx );
                idx += step;
                idx %= N;
            }
            result.push_back(result_idx + 1);
        }
        return result;
    }

    /**
     * Несистематическое декодирование кодом Рида-Соломона. Вектор должен быть скорректирован, 
     * то есть давать нулевой синдром. Фактически, это последний этап декодирования.
     */
    template< int p, int q >
    inline std::vector< int > Decode( const std::vector< int >& v, int R, const gf::GF< p, q >& gf)
    {
        const int N = std::pow( p, q ) - 1;
        std::vector< int > result;
        // a' = v * F', F' - квадратная матрица, обратная матрице F.
        for (int i = 0; i < N; ++i) { // По столбцам матрицы F'.
            const int step = -i; // Шаг степени элементов "альфа".
            int idx = 0;
            int result_idx = -1;
            for (int j = 0; j < N; ++j) { // По строкам i-го столбца.
                const int mult_idx = gf.Mult( v.at( j ) - 1, idx );
                result_idx = gf.Add( mult_idx, result_idx );
                idx += (step + N);
                idx %= N;
            }
            result.push_back(result_idx + 1);
        }
        assert(K > 0);
        const int K = N - R;
        while (result.size() > K) {
            assert(result.back() == 0);
            result.pop_back();
        }
        return result;
    }

    /**
     * Комбинация кода Рида-Соломона (РС) и расширенного кода Хэмминга.
     * Код РС обнаруживает ошибки, и исправляет 1-кратные по полной таблице LUT. В случае
     * невозможности исправить - стирает все текущие символы.
     * Выбран табличный способ из-за быстроты и простоты (по крайней мере для CPU-систем).
     * Чтобы не занимать много памяти, выбран короткий код РС. Плата за это - высоковатая избыточность.
     * Результирующая избыточность (вместе с расширенным кодом Хэмминга) около 44%.
     * Расширенный код Хэмминга работает в режиме восстановления стертых символов. Символом
     * для него является кодовый вектор кода РС. Используется 32 параллельных кода РС.
     */
    struct RsExh {
        static constexpr int p = 2;
        static constexpr int q = 4; // N = p^q - 1 - длина кода Рида-Соломона.
        static constexpr int N = utils::power<int>(p, q) - 1;
        static constexpr int R = 4; // Количество проверочных символов кода Рида-Соломона.
        static constexpr int K = N - R;
        // Кодовое расстояние.
        static constexpr int D = R + 1;
        // Базовая таблица для кода Рида-Соломона. Строится по порождающему полиному.
        gf::GFLUT mLut{ gf::State( p, std::vector< int >{ 1, 0, 0, 1 } ) };
        // Проверка, что полином примитивный.
        const bool mIsGood = mLut.PolyIsGood();
        // Манипулятор полем Галуа GF(p^q).
        gf::GF< p, q > mGf{ mLut };
        // Таблица соответствия синдромов и им соответствующих однократных ошибок (1-ошибок).
        std::unordered_map< std::vector< int >, std::pair<int, int>, gf::KeyHasher2 > mLut_1_errors;
        static constexpr int R2 = 6;  // Количество проверочных символов расширенного кода Хэмминга.
        static constexpr int M2 = 11; // Количество внутренних символов расширенного кода Хэмминга.
        hamming::HammingExtended< R2, M2, int > mHammingCode;

        /**
         * Конструктор. Заполняется таблица для исправления однократных ошибок.
         */
        RsExh()
        {
            assert(mIsGood);
            mLut_1_errors.clear();
            const int N = std::pow( p, q ) - 1;
            mLut_1_errors.reserve( N * (N-1) );
            std::vector< int > c( R ); // Синдром.
            for( int i = 0; i < N; ++i )
            {
                for( int j = 0; j < N - 1; ++j )
                {
                    for( int ii = 0; ii < R; ++ii ) // Строки проверочной матрицы H.
                    {
                        int idx = 0;
                        int result_idx = -1;
                        idx += i * ( ii + 1 );
                        idx %= N;
                        int mult_idx = mGf.Mult( idx + j, 0 ); // ? idx + j, 0
                        result_idx = mGf.Add( mult_idx, result_idx );
                        c[ ii ] = result_idx + 1;
                    }
                    mLut_1_errors[ c ] = std::make_pair(i, j); // Позиция и значение 1-кратной ошибки.
                }
            }
        }
    };
}
