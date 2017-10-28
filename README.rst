Analytic
======
A collection of general purpose analytic solutions for polynomials of degree zero through four in various languages
The roots of a linear, quadratic, cubic and quartic equations can be solved from their coefficients - doing so is generally much faster then numerical methods. This repository aims to be a collection of such methods for use.

Polynomials
-----------
- **degree zero**
  The zero polynomial - a constant of zero
  - a = 0
  
- **degree one**
  The linear equation with constants a, b
  - ax + b = 0

- **degree two**
  The quadratic equation with constants a, b, c
  - ax^2 + bx + c = 0

- **degree three**
  The cubic equation with constants a, b, c, d
  - ax^3 + bx^2 + cx + d = 0

- **degree four**
  The quartic equation with constants a, b, c, d, e
  - ax^4 + bx^3 + cx^2 + dx + e = 0
  
- **higher degrees**
  Higher polynomial degrees do not have analytic solutions

  
Notes
-----
methods for each polynomial should be self-contained. A method of solving for roots of a quartic should not call a method that solves a quadratic or cubic. If it needs to solve a quadratic or cubic, have it do it right there.

The zero polynomial and linear equations are trivial but included for the sake of completion.
