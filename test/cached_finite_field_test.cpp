// Checks that the cached finite field gives the same answers as the original
// polynomial arithmetic it wraps.
//
// F is the cached field under its usual public name. U is the original class,
// kept in the algebra::uncached namespace. The field has 1331 elements, so the
// tests go through every element and every pair of elements instead of sampling.

#include <array>
#include <cstddef>
#include <cstdint>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>
#include <gtest/gtest.h>
#include <algebra/finite_fields/finite_fields.hpp>

using namespace algebra;

namespace {

using F = FiniteField<11, 3>;
using E = F::Element;
using U = uncached::FiniteField<11, 3>;
using UE = U::Element;

static_assert(FieldElement<E>);
static_assert(FieldElement<UE>);
static_assert(E::is_cached);
static_assert(std::is_same_v<E::uncached_type, UE>);
static_assert(std::is_same_v<E::polynomial_type, UE::polynomial_type>);
static_assert(std::is_same_v<E::coefficient_type, UE::coefficient_type>);

// A cached element is a two-byte index that can be copied like an integer.
static_assert(sizeof(E) == sizeof(std::uint16_t));
static_assert(std::is_trivially_copyable_v<E>);

// 4096 elements is the largest field that gets tables.
static_assert(FiniteField<2, 12>::Element::is_cached);
static_assert(!FiniteField<2, 13>::Element::is_cached);
static_assert(!FiniteField<2, 64>::Element::is_cached);
static_assert(FieldElement<FiniteField<2, 13>::Element>);

template<class Field>
std::vector<typename Field::Element> every_element() {
    std::vector<typename Field::Element> result;
    const std::uint64_t count = Field::cardinality();
    result.reserve(static_cast<std::size_t>(count));
    for (std::uint64_t index = 0; index < count; ++index) {
        result.push_back(Field::Element::from_index(index));
    }
    return result;
}

template<class Value>
std::string text(const Value& value) {
    std::ostringstream out;
    out << value;
    return out.str();
}

// Runs one binary operation on every ordered pair of elements, twice. The
// first pass fills the table from the original arithmetic and the second pass
// reads every answer back from the table.
template<class Cached, class Original>
void expect_same_for_every_pair(Cached cached_operation, Original original_operation) {
    const auto cached = every_element<F>();
    const auto original = every_element<U>();
    for (int pass = 0; pass < 2; ++pass) {
        std::uint64_t mismatches = 0;
        for (std::size_t a = 0; a < cached.size(); ++a) {
            for (std::size_t b = 0; b < cached.size(); ++b) {
                const E actual = cached_operation(cached[a], cached[b]);
                const UE expected = original_operation(original[a], original[b]);
                if (!(actual.representative() == expected.representative())) {
                    ++mismatches;
                }
            }
        }
        EXPECT_EQ(mismatches, 0u) << "pass " << pass;
    }
}

TEST(CachedFiniteField, DescribesTheSameFieldAsTheOriginal) {
    EXPECT_EQ(&F::modulus(), &U::modulus());
    EXPECT_EQ(F::cardinality(), U::cardinality());
    EXPECT_EQ(F::characteristic, U::characteristic);
    EXPECT_EQ(F::extension_degree, U::extension_degree);
    EXPECT_EQ(F::cache_statistics().order, 1331u);
    EXPECT_EQ(F::cache_statistics().table_bytes, 2u * (2u * 1331u * 1331u + 2u * 1331u));
}

TEST(CachedFiniteField, NumbersEveryElementLikeTheOriginal) {
    const auto cached = every_element<F>();
    const auto original = every_element<U>();
    ASSERT_EQ(cached.size(), 1331u);
    ASSERT_EQ(original.size(), 1331u);
    for (std::size_t i = 0; i < cached.size(); ++i) {
        ASSERT_EQ(cached[i].representative(), original[i].representative()) << i;
        ASSERT_EQ(cached[i].coefficients(), original[i].coefficients()) << i;
        ASSERT_EQ(cached[i].is_zero(), original[i].is_zero()) << i;
        ASSERT_EQ(cached[i].is_one(), original[i].is_one()) << i;
        ASSERT_EQ(text(cached[i]), text(original[i])) << i;
        for (std::size_t j = 0; j < cached.size(); ++j) {
            ASSERT_EQ(cached[i] == cached[j], original[i] == original[j]) << i << ' ' << j;
        }
    }
    EXPECT_THROW(E::from_index(1331), std::out_of_range);
    EXPECT_THROW(UE::from_index(1331), std::out_of_range);
    EXPECT_THROW(E::one().coefficient(3), std::out_of_range);
}

TEST(CachedFiniteField, AddsEveryPairLikeTheOriginal) {
    expect_same_for_every_pair([](const E& a, const E& b) { return a.add(b); },
                               [](const UE& a, const UE& b) { return a.add(b); });
    EXPECT_EQ(F::cache_statistics().sums, 1331u * 1331u);
}

TEST(CachedFiniteField, MultipliesEveryPairLikeTheOriginal) {
    expect_same_for_every_pair([](const E& a, const E& b) { return a.multiply(b); },
                               [](const UE& a, const UE& b) { return a.multiply(b); });
    EXPECT_EQ(F::cache_statistics().products, 1331u * 1331u);
}

TEST(CachedFiniteField, NegatesAndInvertsEveryElementLikeTheOriginal) {
    const auto cached = every_element<F>();
    const auto original = every_element<U>();
    for (int pass = 0; pass < 2; ++pass) {
        for (std::size_t i = 0; i < cached.size(); ++i) {
            ASSERT_EQ(cached[i].negate().representative(), original[i].negate().representative())
                << i;
            if (i == 0) {
                EXPECT_THROW(cached[i].reciprocal(), std::domain_error);
                EXPECT_THROW(original[i].reciprocal(), std::domain_error);
            } else {
                ASSERT_EQ(cached[i].reciprocal().representative(),
                          original[i].reciprocal().representative()) << i;
            }
        }
    }
    EXPECT_EQ(F::cache_statistics().negatives, 1331u);
    EXPECT_EQ(F::cache_statistics().reciprocals, 1330u);
}

TEST(CachedFiniteField, OperatorsAgreeWithTheOriginal) {
    const auto cached = every_element<F>();
    const auto original = every_element<U>();
    std::mt19937 generator(2026u);
    for (int trial = 0; trial < 20000; ++trial) {
        const std::size_t a = generator() % cached.size();
        const std::size_t b = generator() % cached.size();
        ASSERT_EQ((cached[a] + cached[b]).representative(),
                  (original[a] + original[b]).representative());
        ASSERT_EQ((cached[a] - cached[b]).representative(),
                  (original[a] - original[b]).representative());
        ASSERT_EQ((-cached[a]).representative(), (-original[a]).representative());
        ASSERT_EQ((cached[a] * cached[b]).representative(),
                  (original[a] * original[b]).representative());
        if (b != 0) {
            ASSERT_EQ((cached[a] / cached[b]).representative(),
                      (original[a] / original[b]).representative());
        }
        E accumulated = cached[a];
        UE expected = original[a];
        accumulated += cached[b];
        expected += original[b];
        accumulated *= cached[a];
        expected *= original[a];
        accumulated -= cached[b];
        expected -= original[b];
        ASSERT_EQ(accumulated.representative(), expected.representative());
    }
}

TEST(CachedFiniteField, BuildsElementsLikeTheOriginal) {
    EXPECT_EQ(E{}.representative(), UE{}.representative());
    EXPECT_EQ(F::zero().representative(), U::zero().representative());
    EXPECT_EQ(F::one().representative(), U::one().representative());
    EXPECT_EQ(F::alpha().representative(), U::alpha().representative());
    for (std::int64_t value = -40; value <= 40; ++value) {
        EXPECT_EQ(E::from_integer(value).representative(),
                  UE::from_integer(value).representative());
        EXPECT_EQ(F::element(value).representative(), U::element(value).representative());
    }
    // Coefficient lists longer than the extension degree have to be reduced by the modulus.
    EXPECT_EQ(F::element({1, 2, 3, 4, 5}).representative(),
              U::element({1, 2, 3, 4, 5}).representative());
    EXPECT_EQ(F::element({-1, 25, -7}).representative(), U::element({-1, 25, -7}).representative());
    EXPECT_EQ(E::from_coefficients({0, 0, 0, 1}).representative(),
              UE::from_coefficients({0, 0, 0, 1}).representative());
    std::mt19937 generator(7u);
    for (int trial = 0; trial < 2000; ++trial) {
        F::PolynomialType polynomial;
        const std::size_t length = 1 + generator() % 9;
        for (std::size_t degree = 0; degree < length; ++degree) {
            const auto coefficient = static_cast<std::int64_t>(generator() % 11);
            polynomial.set_coefficient(degree, E::coefficient_type::from_integer(coefficient));
        }
        ASSERT_EQ(E::from_polynomial(polynomial).representative(),
                  UE::from_polynomial(polynomial).representative());
    }
}

// The same product and elimination the assignment program runs, written once
// for any element type so it can be run on both implementations.
template<class T>
std::vector<T> random_matrix(std::size_t n, std::uint32_t seed) {
    std::mt19937 generator(seed);
    std::vector<T> matrix;
    matrix.reserve(n * n);
    for (std::size_t i = 0; i < n * n; ++i) {
        matrix.push_back(T::from_index(generator() % 1331));
    }
    return matrix;
}

template<class T>
std::vector<T> product(const std::vector<T>& a, const std::vector<T>& b, std::size_t n) {
    std::vector<T> result(n * n, T::zero());
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < n; ++j) {
            T sum = T::zero();
            for (std::size_t k = 0; k < n; ++k) {
                sum = sum.add(a[i * n + k].multiply(b[k * n + j]));
            }
            result[i * n + j] = sum;
        }
    }
    return result;
}

template<class T>
T determinant(std::vector<T> matrix, std::size_t n) {
    T result = T::one();
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot_row = k;
        while (pivot_row < n && matrix[pivot_row * n + k].is_zero()) {
            ++pivot_row;
        }
        if (pivot_row == n) {
            return T::zero();
        }
        if (pivot_row != k) {
            for (std::size_t j = 0; j < n; ++j) {
                std::swap(matrix[pivot_row * n + j], matrix[k * n + j]);
            }
            result = result.negate();
        }
        const T pivot = matrix[k * n + k];
        result = result.multiply(pivot);
        T inverse = T::one();
        if (!pivot.is_one()) {
            inverse = pivot.reciprocal();
        }
        for (std::size_t i = k + 1; i < n; ++i) {
            if (matrix[i * n + k].is_zero()) {
                continue;
            }
            const T multiplier = matrix[i * n + k].multiply(inverse);
            matrix[i * n + k] = T::zero();
            for (std::size_t j = k + 1; j < n; ++j) {
                const T correction = multiplier.multiply(matrix[k * n + j]);
                matrix[i * n + j] = matrix[i * n + j].add(correction.negate());
            }
        }
    }
    return result;
}

TEST(CachedFiniteField, ComputesTheSameDeterminantsAsTheOriginal) {
    constexpr std::size_t n = 40;
    for (const std::uint32_t seed : {123456u, 7u, 99u}) {
        const auto a = random_matrix<E>(n, seed);
        const auto b = random_matrix<E>(n, seed + 1);
        const auto original_a = random_matrix<UE>(n, seed);
        const auto original_b = random_matrix<UE>(n, seed + 1);

        const E det_a = determinant(a, n);
        const E det_b = determinant(b, n);
        const E det_c = determinant(product(a, b, n), n);
        const UE original_det_a = determinant(original_a, n);
        const UE original_det_b = determinant(original_b, n);
        const UE original_det_c = determinant(product(original_a, original_b, n), n);

        EXPECT_EQ(det_a.representative(), original_det_a.representative());
        EXPECT_EQ(det_b.representative(), original_det_b.representative());
        EXPECT_EQ(det_c.representative(), original_det_c.representative());
        EXPECT_EQ(det_c, det_a.multiply(det_b));
    }
}

TEST(CachedFiniteField, FillsTheTablesCorrectlyFromSeveralThreads) {
    // A field no other test touches, so every thread starts on empty tables.
    using Small = FiniteField<7, 3>;
    using SmallOriginal = uncached::FiniteField<7, 3>;
    Small::initialize();
    const auto cached = every_element<Small>();
    const auto original = every_element<SmallOriginal>();
    const std::size_t count = cached.size();
    ASSERT_EQ(count, 343u);

    std::vector<std::uint64_t> mismatches(4, 0);
    std::vector<std::thread> threads;
    for (std::size_t id = 0; id < mismatches.size(); ++id) {
        threads.emplace_back([&, id] {
            for (std::size_t a = 0; a < count; ++a) {
                for (std::size_t b = 0; b < count; ++b) {
                    const auto sum = cached[a].add(cached[b]);
                    const auto times = cached[a].multiply(cached[b]);
                    if (!(sum.representative() == original[a].add(original[b]).representative()) ||
                        !(times.representative() ==
                          original[a].multiply(original[b]).representative())) {
                        ++mismatches[id];
                    }
                }
            }
        });
    }
    for (auto& thread : threads) {
        thread.join();
    }
    for (const auto count_for_thread : mismatches) {
        EXPECT_EQ(count_for_thread, 0u);
    }
    EXPECT_EQ(Small::cache_statistics().sums, 343u * 343u);
    EXPECT_EQ(Small::cache_statistics().products, 343u * 343u);
}

TEST(CachedFiniteField, ForwardsFieldsThatAreTooLargeToTabulate) {
    using Big = FiniteField<2, 13>;
    using BigOriginal = uncached::FiniteField<2, 13>;
    EXPECT_EQ(Big::cardinality(), 8192u);
    EXPECT_EQ(Big::modulus(), BigOriginal::modulus());
    EXPECT_EQ(Big::cache_statistics().order, 0u);
    EXPECT_EQ(Big::cache_statistics().table_bytes, 0u);
    std::mt19937 generator(13u);
    for (int trial = 0; trial < 500; ++trial) {
        const std::uint64_t a = generator() % 8192;
        const std::uint64_t b = 1 + generator() % 8191;
        const auto x = Big::Element::from_index(a);
        const auto y = Big::Element::from_index(b);
        const auto original_x = BigOriginal::Element::from_index(a);
        const auto original_y = BigOriginal::Element::from_index(b);
        ASSERT_EQ((x + y).representative(), (original_x + original_y).representative());
        ASSERT_EQ((x * y).representative(), (original_x * original_y).representative());
        ASSERT_EQ((-x).representative(), (-original_x).representative());
        ASSERT_EQ(y.reciprocal().representative(), original_y.reciprocal().representative());
        ASSERT_TRUE((y * y.reciprocal()).is_one());
    }
    EXPECT_TRUE(Big::zero().is_zero());
    EXPECT_TRUE(Big::one().is_one());
    EXPECT_THROW(Big::zero().reciprocal(), std::domain_error);
}

}
