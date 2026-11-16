# Algebra finite fields (cached)

A drop-in replacement for [algebra-finite-fields](https://github.com/nathanjrussell/algebra-finite-fields) that caches the results of finite field arithmetic. It is based on upstream commit `74b3d56` and keeps the same CMake target (`algebra::finite_fields`), header paths and type names, so a project switches to it by changing one dependency.

## Using it

In the project that depends on `algebra-finite-fields`, point the existing `algebra_finite_fields` entry at this repository. Nothing else changes.

```cmake
FetchContent_Declare(
    algebra_finite_fields
    GIT_REPOSITORY <URL of this repository>
    GIT_TAG <commit to use>
)
```

## How the cache works

`FiniteField<11, 3>` has 11³ = 1331 elements. The original class stores each element as a polynomial and repeats the polynomial multiplication and reduction on every operation. The cached class stores each element as a two-byte index, its position in the `from_index` numbering, and keeps one table per operation.

| Table | Slots | A slot holds |
| --- | --- | --- |
| sums | q × q | the index of a + b |
| products | q × q | the index of a × b |
| negatives | q | the index of −a |
| reciprocals | q | the index of 1 / a |

A slot is empty until that operation is first requested. On that first request, a cache miss, the result is computed by the original class and stored. Every later request for the same operands, a cache hit, is one table read. For q = 1331 the four tables take 7.1 MB.

The original classes are still in the library as `algebra::uncached::FiniteField` and `algebra::uncached::FiniteFieldElement` in `uncached_finite_field.hpp`. The cached field uses the modulus the original class selects, so both describe the same field.

## What changed from upstream

| File | Change |
| --- | --- |
| `include/algebra/finite_fields/uncached_finite_field.hpp` | The upstream `finite_field.hpp`, renamed. The only edit is the namespace, `algebra` to `algebra::uncached`. |
| `include/algebra/finite_fields/finite_field.hpp` | New. The cached `FiniteField` and `FiniteFieldElement`. |
| `test/cached_finite_field_test.cpp` | New. Compares the cached class with the original. |
| `test/CMakeLists.txt` | Builds the new test file. |

`modulus.hpp`, `prime_field.hpp`, `enumeration.hpp`, the upstream test file and the other build files are unchanged.

## Tests

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

- `test/finite_field_test.cpp` is the upstream test file, unmodified. It now runs against the cached class.
- `test/cached_finite_field_test.cpp` compares the cached class with the original for every element and every ordered pair of elements, which is 1331² = 1,771,561 pairs each for add and multiply. Each pair is checked twice so that both the miss and the hit are covered. It also runs the same matrix product and elimination as the assignment program on both implementations and checks that the determinants match.

## Limits

- Fields with more than 4096 elements are not cached. Their elements hold the original element and forward every call to it. The limit is `ALGEBRA_FINITE_FIELDS_CACHE_LIMIT`, and the two pair tables take 4q² bytes.
- The tables can be read and filled from several threads at once. Slots are atomic and a slot only ever receives one value.
- `FiniteField::cache_statistics()`, `FiniteFieldElement::is_cached` and `FiniteFieldElement::uncached_type` are additions to the upstream interface. Nothing was removed.
