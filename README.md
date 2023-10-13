# HugeNumbers.jl

A package for working with **huge or tiny numbers** beyond the dynamic range of
ordinary floating point numbers.

Exports the number type `Huge{T} <: Real` which is a lot like a floating point number,
except it has an exponentially larger range.

You can use this for calculating extremely small probabilities, which would normally
underflow a float, or extremely large quantities which would normally overflow.

## Usage

To convert a number `x` to a `Huge` number:
- Do `convert(Huge, x)` or `Huge(x)`.
- If you already know `ix = hlog(x)` then `hexp(Huge, ix)` is faster.
- If you already know `lx = log(x)` then `exp(Huge, lx)` is faster.

Convert a number `x` to a `Huge` number with `convert(Huge, x)` or `Huge(x)`. If you
happen to already know `ix = hlog(x)` then `hexp(Huge, ix)` is faster.

The following operations are implemented and accurate:
- **Arithmetic:** `+`, `-`, `*`, `/`, `^`, `inv`.
- **Ordering:** `==`, `<`, `cmp`, `isequal`, `isless`, `sign`, `signbit`, `abs`.
- **Logarithm:** `log`, `log2`, `log10`].
- **Conversion:** `float`, `widen`, `big`.
- **Special values:** `zero`, `one`, `typemin`, `typemax`.
- **Predicates:** `iszero`, `isone`, `isinf`, `isfinite`, `isnan`.
- **IO:** `show`, `read`, `write`.
- **Misc:** `nextfloat`, `prevfloat`, `hash`.
- **Note:** Any functions not mentioned here might be inaccurate.

### Interoperability with other packages

It is natural to use this package in conjunction with other packages which return
logarithms. The general pattern is that you can use `exp(Huge, logfunc(args...))`
instead of `func(args...)` to get the answer as a logarithmic number. Here are some
possibilities for `func`:

- [StatsFuns.jl](https://github.com/JuliaStats/StatsFuns.jl):
  `normpdf`, `normcdf`, `normccdf`, plus equivalents for other distributions.
- [Distributions.jl](https://github.com/JuliaStats/Distributions.jl):
  `pdf`, `cdf`, `ccdf`.
- [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl):
  `gamma`, `factorial`, `beta`, `erfc`, `erfcx`.

## Relationship to [LogarithmicNumbers](https://github.com/cjdoris/LogarithmicNumbers.jl)

Whereas a logarithmic number stores `log(x)`, we store a different quantity `hlog(x)`.
This is the inverse of `hexp(x)`, which has the following nice properties:
- It is an increasing, continuous, differentiable, invertible function on the real numbers.
- It grows exponentially for large `x` and shrinks exponentially for small `x`.
- `hexp(x) = x` for `x` in `-Inf`, `-1`, `0`, `1` and `Inf`.
- `-hexp(x) = hexp(-x)`.
- `1/hexp(x) = hexp(1/x)`.

The crucial difference between `log(x)` and `hlog(x)` is that the former is only valid
for positive numbers. This means that a logrithmic number requires an extra sign bit if you
need to store a sign, whereas a `Huge{T}` can be negative provided `T` can be negative.

The other properties mean that many arthmetic and comparison operations are trivial to
implement - often much simpler than with logarithmic numbers.

`hexp(x)` is defined as:
- `exp(x - 1)` if `x ≥ 1`
- `exp(1 - 1/x)` if `x ≥ 0`
- `-exp(1 + 1/x)` if `x ≥ -1`
- `-exp(-x - 1)` otherwise
