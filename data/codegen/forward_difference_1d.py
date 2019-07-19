#!/usr/bin/env python3

import numpy as np
import scipy.special
from sympy import *
from enum import Enum
import sys
import re

x,h,f = symbols('x,h,f')

MAX_DEGREE = 12

monomials = []
coeffs = []

for d in range(MAX_DEGREE+1):
  monomials.append(x**d)
  coeffs.append(symbols('a{}'.format(d)))
  coeffs[d].name = 'poly1_.getCoefficient({})'.format(d)

def replacePow(string):

  def process_match(m):
    symbol = m[1]
    exponent = int(m[2])

    s = symbol
    for k in range(1, exponent):
      s += '*' + symbol

    return s


  p = re.compile("pow\((\w), (\d+)\)")
  string = p.sub(process_match, string)

  return string

def generateStartX(degree):
  f = 0

  for k in range(degree+1):
    f += coeffs[k] * monomials[k]

  delta = horner(f.subs(x, x + h) - f)

  print('    deltas[{}] = {};'.format(0, replacePow(ccode(f))))

  if (degree >= 1):
    print('    deltas[{}] = {};'.format(1, replacePow(ccode(delta))))

  for k in range(2, degree+1):  
    delta = horner(delta.subs(x, x+h) - delta)
    print('     deltas[{}] = {};'.format(k, replacePow(ccode(delta))))

print('real_t setStartX(real_t x, real_t h) {')
print('  static_assert(Degree <= {}, "Polynomial2DEvaluator not implemented for degree larger than {}");'.format(MAX_DEGREE, MAX_DEGREE))

for degree in range(MAX_DEGREE+1):
  print('  if (Degree == {}) {{'.format(degree))
  generateStartX(degree)
  print('  }')

print('  return deltas[0];')
print('}\n')

print('real_t incrementEval() {')

for k in range(1, MAX_DEGREE+1):  
  print('  if (Degree >= {}) {{'.format(k))
  print('     deltas[{}] += deltas[{}];'.format(k-1, k))
  print('  }')

print('  return deltas[0];')
print('}')
