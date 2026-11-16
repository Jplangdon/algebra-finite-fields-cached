#pragma once

// Cached finite field.
//
// This header keeps the public names of the original library
// (algebra::FiniteField and algebra::FiniteFieldElement) but puts a cache in
// front of the arithmetic. The original classes still exist, unchanged apart
// from their namespace, as algebra::uncached::FiniteField and
// algebra::uncached::FiniteFieldElement in uncached_finite_field.hpp.
//
// How it works
//   A field with q = P^K elements only has q possible values, so an element is
//   stored as a small integer, its index, instead of as a polynomial. The
//   index uses the same numbering as from_index: index = c0 + c1*P + c2*P^2 +
//   ... where c0, c1, ... are the polynomial coefficients.
//
//   Each operation has a table with one slot per operand (negate, reciprocal)
//   or per pair of operands (add, multiply). A slot starts out empty. The
//   first time an operation is asked for (a cache miss) the answer is computed
//   by the original polynomial arithmetic and written into the slot. Every
//   later request for the same operands (a cache hit) is a single table read.
//
// Fields that are too large to tabulate (more than
// ALGEBRA_FINITE_FIELDS_CACHE_LIMIT elements) skip the cache: their elements
// hold the original element and forward every call to it.

#include <array>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <ostream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>
#include <algebra/finite_fields/modulus.hpp>
#include <algebra/finite_fields/uncached_finite_field.hpp>

// Largest field, counted in elements, that gets tables. The two pair tables
// together take 4 * q^2 bytes: about 7 MB for q = 1331 and 67 MB for q = 4096.
#ifndef ALGEBRA_FINITE_FIELDS_CACHE_LIMIT
#define ALGEBRA_FINITE_FIELDS_CACHE_LIMIT 4096
#endif

// Marks the functions that run on a cache miss so the compiler keeps them out
// of the callers' loops.
#if defined(__GNUC__) || defined(__clang__)
#define ALGEBRA_FINITE_FIELDS_CACHE_MISS __attribute__((noinline, cold))
#else
#define ALGEBRA_FINITE_FIELDS_CACHE_MISS
#endif

namespace algebra {

// Number of results currently held in each table. Returned by
// FiniteField::cache_statistics().
struct FiniteFieldCacheStatistics {
    std::uint64_t order = 0;         // elements in the field, or 0 if the field is not cached
    std::uint64_t sums = 0;          // filled slots of the addition table, at most order^2
    std::uint64_t products = 0;      // filled slots of the multiplication table, at most order^2
    std::uint64_t negatives = 0;     // filled slots of the negation table, at most order
    std::uint64_t reciprocals = 0;   // filled slots of the reciprocal table, at most order - 1
    std::uint64_t table_bytes = 0;   // size of the four tables
};

namespace detail {

inline constexpr std::uint64_t finite_field_cache_limit = ALGEBRA_FINITE_FIELDS_CACHE_LIMIT;
static_assert(finite_field_cache_limit < 65535,
              "ALGEBRA_FINITE_FIELDS_CACHE_LIMIT must leave room for a 16-bit element index");

// base^exponent if that is at most limit, otherwise 0. Stops before it can overflow.
constexpr std::uint64_t bounded_power(std::uint64_t base, std::size_t exponent,
                                      std::uint64_t limit) {
    std::uint64_t result = 1;
    for (std::size_t i = 0; i < exponent; ++i) {
        if (result > limit / base) {
            return 0;
        }
        result *= base;
    }
    return result;
}

// P^K for fields small enough to cache, 0 for every other field.
template<std::uint32_t P, std::size_t K>
inline constexpr std::uint64_t cached_field_order =
    bounded_power(P, K, finite_field_cache_limit);

// The tables for one field.
template<std::uint32_t P, std::size_t K, class ModulusPolicy>
class FiniteFieldCache {
public:
    using element_type = uncached::FiniteFieldElement<P, K, ModulusPolicy>;
    using index_type = std::uint16_t;
    static constexpr std::size_t order = static_cast<std::size_t>(cached_field_order<P, K>);

private:
    static_assert(order > 0 && order < std::numeric_limits<index_type>::max());

    // A slot holds 0 until its result has been computed. After that it holds
    // the index of the result plus one.
    //
    // Slots are atomic so that two threads may use the field at the same time.
    // If both miss on the same slot they compute the same value and store the
    // same value, so no lock is needed. A relaxed load is an ordinary read.
    using slot_type = std::atomic<index_type>;

    // The tables have static storage, which starts out as zeros. Every slot is
    // therefore empty at program start, there is nothing to set up, and the
    // operating system only backs the pages that are actually written.
    static inline slot_type sums_[order * order]{};       // slot a * order + b is for a + b
    static inline slot_type products_[order * order]{};   // slot a * order + b is for a * b
    static inline slot_type negatives_[order]{};          // slot a is for -a
    static inline slot_type reciprocals_[order]{};        // slot a is for 1 / a

    // The original elements, used to compute a result on a cache miss.
    struct Originals {
        std::vector<element_type> elements;   // elements[i] is the original element with index i
        index_type one = 0;

        Originals() {
            elements.reserve(order);
            for (std::size_t index = 0; index < order; ++index) {
                elements.push_back(element_type::from_index(index));
                if (index_of(elements.back()) != index) {
                    throw std::logic_error(
                        "Finite field indices do not match the element numbering");
                }
            }
            one = index_of(element_type::one());
        }
    };

    static const Originals& originals() {
        static const Originals instance;
        return instance;
    }

    // Cache misses. Each one runs the original arithmetic once and remembers
    // the result. They are kept out of line so that the table reads in add(),
    // multiply(), negate() and reciprocal() stay small enough to inline.
    static index_type store(slot_type& slot, const element_type& result) {
        const index_type index = index_of(result);
        slot.store(static_cast<index_type>(index + 1), std::memory_order_relaxed);
        return index;
    }

    ALGEBRA_FINITE_FIELDS_CACHE_MISS static index_type add_miss(index_type a, index_type b) {
        const auto& elements = originals().elements;
        return store(sums_[std::size_t{a} * order + b], elements[a].add(elements[b]));
    }

    ALGEBRA_FINITE_FIELDS_CACHE_MISS static index_type multiply_miss(index_type a, index_type b) {
        const auto& elements = originals().elements;
        return store(products_[std::size_t{a} * order + b], elements[a].multiply(elements[b]));
    }

    ALGEBRA_FINITE_FIELDS_CACHE_MISS static index_type negate_miss(index_type a) {
        return store(negatives_[a], originals().elements[a].negate());
    }

    // The original reciprocal throws for zero, so nothing is stored for it.
    ALGEBRA_FINITE_FIELDS_CACHE_MISS static index_type reciprocal_miss(index_type a) {
        return store(reciprocals_[a], originals().elements[a].reciprocal());
    }

    template<std::size_t Size>
    static std::uint64_t filled(const slot_type (&table)[Size]) {
        std::uint64_t count = 0;
        for (const slot_type& slot : table) {
            if (slot.load(std::memory_order_relaxed) != 0) {
                ++count;
            }
        }
        return count;
    }

public:
    // Builds the list of original elements, which also checks the modulus.
    static void initialize() {
        static_cast<void>(originals());
    }

    // Index of an original element: its coefficients read as a base-P number.
    static index_type index_of(const element_type& value) {
        std::size_t index = 0;
        for (std::size_t i = K; i > 0; --i) {
            index = index * P + value.coefficient(i - 1).representative();
        }
        return static_cast<index_type>(index);
    }

    static const element_type& element(index_type index) {
        return originals().elements[index];
    }

    static index_type one() {
        return originals().one;
    }

    static index_type add(index_type a, index_type b) {
        const index_type slot = sums_[std::size_t{a} * order + b].load(std::memory_order_relaxed);
        if (slot != 0) [[likely]] {
            return static_cast<index_type>(slot - 1);
        }
        return add_miss(a, b);
    }

    static index_type multiply(index_type a, index_type b) {
        const index_type slot =
            products_[std::size_t{a} * order + b].load(std::memory_order_relaxed);
        if (slot != 0) [[likely]] {
            return static_cast<index_type>(slot - 1);
        }
        return multiply_miss(a, b);
    }

    static index_type negate(index_type a) {
        const index_type slot = negatives_[a].load(std::memory_order_relaxed);
        if (slot != 0) [[likely]] {
            return static_cast<index_type>(slot - 1);
        }
        return negate_miss(a);
    }

    static index_type reciprocal(index_type a) {
        const index_type slot = reciprocals_[a].load(std::memory_order_relaxed);
        if (slot != 0) [[likely]] {
            return static_cast<index_type>(slot - 1);
        }
        return reciprocal_miss(a);
    }

    static FiniteFieldCacheStatistics statistics() {
        FiniteFieldCacheStatistics result;
        result.order = order;
        result.sums = filled(sums_);
        result.products = filled(products_);
        result.negatives = filled(negatives_);
        result.reciprocals = filled(reciprocals_);
        result.table_bytes =
            sizeof(sums_) + sizeof(products_) + sizeof(negatives_) + sizeof(reciprocals_);
        return result;
    }
};

}

template<std::uint32_t P, std::size_t K = 1, class ModulusPolicy = FirstEnumeratedIrreducible>
class FiniteField;

template<std::uint32_t P, std::size_t K, class ModulusPolicy>
class FiniteFieldElement : public AdditiveOperators<FiniteFieldElement<P, K, ModulusPolicy>>,
                           public MultiplicativeOperators<FiniteFieldElement<P, K, ModulusPolicy>>,
                           public DivisionOperators<FiniteFieldElement<P, K, ModulusPolicy>> {
    static_assert(is_prime(P));
    static_assert(K > 0);

public:
    using parent_type = FiniteField<P, K, ModulusPolicy>;
    using algebra = FieldDescriptor;
    using coefficient_type = PrimeFieldElement<P>;
    using polynomial_type = Polynomial<coefficient_type>;

    // The original element type that does the arithmetic on a cache miss.
    using uncached_type = uncached::FiniteFieldElement<P, K, ModulusPolicy>;

    // True when this field is small enough to use the tables.
    static constexpr bool is_cached = detail::cached_field_order<P, K> != 0;

private:
    using cache = detail::FiniteFieldCache<P, K, ModulusPolicy>;

    // A cached field stores the element's index (two bytes, zero is index 0).
    // Any other field stores the original element itself.
    using storage_type = std::conditional_t<is_cached, std::uint16_t, uncached_type>;
    storage_type value_{};
    struct Stored {};

    FiniteFieldElement(storage_type value, Stored) : value_(std::move(value)) {}

    static FiniteFieldElement wrap(const uncached_type& value) {
        if constexpr (is_cached) {
            return FiniteFieldElement(cache::index_of(value), Stored{});
        } else {
            return FiniteFieldElement(value, Stored{});
        }
    }

    const uncached_type& unwrap() const {
        if constexpr (is_cached) {
            return cache::element(value_);
        } else {
            return value_;
        }
    }

public:
    FiniteFieldElement() {
        static_cast<void>(parent_type::modulus());
    }

    static FiniteFieldElement from_polynomial(const polynomial_type& value) {
        return wrap(uncached_type::from_polynomial(value));
    }

    static FiniteFieldElement from_coefficients(std::initializer_list<std::int64_t> values) {
        return wrap(uncached_type::from_coefficients(values));
    }

    static FiniteFieldElement from_integer(std::int64_t value) {
        return wrap(uncached_type::from_integer(value));
    }

    static FiniteFieldElement from_index(std::uint64_t index) {
        return wrap(uncached_type::from_index(index));
    }

    static FiniteFieldElement zero() {
        return FiniteFieldElement{};
    }

    static FiniteFieldElement one() {
        if constexpr (is_cached) {
            return FiniteFieldElement(cache::one(), Stored{});
        } else {
            return wrap(uncached_type::one());
        }
    }

    FiniteFieldElement add(const FiniteFieldElement& other) const {
        if constexpr (is_cached) {
            return FiniteFieldElement(cache::add(value_, other.value_), Stored{});
        } else {
            return FiniteFieldElement(value_.add(other.value_), Stored{});
        }
    }

    FiniteFieldElement negate() const {
        if constexpr (is_cached) {
            return FiniteFieldElement(cache::negate(value_), Stored{});
        } else {
            return FiniteFieldElement(value_.negate(), Stored{});
        }
    }

    FiniteFieldElement multiply(const FiniteFieldElement& other) const {
        if constexpr (is_cached) {
            return FiniteFieldElement(cache::multiply(value_, other.value_), Stored{});
        } else {
            return FiniteFieldElement(value_.multiply(other.value_), Stored{});
        }
    }

    FiniteFieldElement reciprocal() const {
        if constexpr (is_cached) {
            return FiniteFieldElement(cache::reciprocal(value_), Stored{});
        } else {
            return FiniteFieldElement(value_.reciprocal(), Stored{});
        }
    }

    bool is_zero() const {
        if constexpr (is_cached) {
            return value_ == 0;
        } else {
            return value_.is_zero();
        }
    }

    bool is_one() const {
        if constexpr (is_cached) {
            return value_ == cache::one();
        } else {
            return value_.is_one();
        }
    }

    coefficient_type coefficient(std::size_t index) const {
        return unwrap().coefficient(index);
    }

    std::array<coefficient_type, K> coefficients() const {
        return unwrap().coefficients();
    }

    const polynomial_type& representative() const {
        return unwrap().representative();
    }

    const parent_type& parent() const {
        static const parent_type instance;
        return instance;
    }

    bool operator==(const FiniteFieldElement& other) const {
        return value_ == other.value_;
    }

    friend std::ostream& operator<<(std::ostream& out, const FiniteFieldElement& value) {
        return out << value.unwrap();
    }
};

template<std::uint32_t P, std::size_t K, class ModulusPolicy>
class FiniteField : public FieldStructure<FiniteFieldElement<P, K, ModulusPolicy>> {
    static_assert(is_prime(P));
    static_assert(K > 0);

    using uncached_field = uncached::FiniteField<P, K, ModulusPolicy>;

public:
    using Element = FiniteFieldElement<P, K, ModulusPolicy>;
    using Base = PrimeField<Cyclic<P>>;
    using PolynomialType = typename Element::polynomial_type;
    static constexpr std::uint32_t characteristic = P;
    static constexpr std::size_t extension_degree = K;

    FiniteField() {
        initialize();
    }

    // The modulus is the one the original class selects, so both
    // implementations describe the same field.
    static const PolynomialType& modulus() {
        return uncached_field::modulus();
    }

    static void initialize() {
        uncached_field::initialize();
        if constexpr (Element::is_cached) {
            detail::FiniteFieldCache<P, K, ModulusPolicy>::initialize();
        }
    }

    static std::uint64_t cardinality() {
        return uncached_field::cardinality();
    }

    static Element element(std::initializer_list<std::int64_t> values) {
        return Element::from_coefficients(values);
    }

    static Element element(std::int64_t value) {
        return Element::from_integer(value);
    }

    static Element zero() {
        return Element::zero();
    }

    static Element one() {
        return Element::one();
    }

    static Element alpha() {
        return Element::from_polynomial(PolynomialType::variable());
    }

    // How much of each table has been filled so far. All zeros for a field
    // that is not cached.
    static FiniteFieldCacheStatistics cache_statistics() {
        if constexpr (Element::is_cached) {
            return detail::FiniteFieldCache<P, K, ModulusPolicy>::statistics();
        } else {
            return FiniteFieldCacheStatistics{};
        }
    }
};

}
