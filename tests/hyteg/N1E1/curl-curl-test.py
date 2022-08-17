import math
import sympy as sp
from sympy.abc import x, y, z

def curl(v):
  return sp.Matrix([ v[2].diff(y) - v[1].diff(z)
                   , v[0].diff(z) - v[2].diff(x)
                   , v[1].diff(x) - v[0].diff(y)
                   ])

def basis(verts):
  canonical_basis = sp.Matrix([ [ 1,  0,  0]
                              , [ 0,  1,  0]
                              , [ 0,  0,  1]
                              , [ 0, -z,  y]
                              , [ z,  0, -x]
                              , [-y,  x,  0]
                              ]
                             ).T

  edges_K_ref = [ (verts[2], verts[3])
                , (verts[1], verts[3])
                , (verts[1], verts[2])
                , (verts[0], verts[3])
                , (verts[0], verts[2])
                , (verts[0], verts[1])
                ]

  dir_K_ref       = list(map(lambda v:  v[1] - v[0]     , edges_K_ref))
  midpoints_K_ref = list(map(lambda v: (v[1] + v[0]) / 2, edges_K_ref))

  A = sp.zeros(6, 6)
  for i in  range(6):
    for j in range(6):
      A[i, j] = canonical_basis.col(j).dot(dir_K_ref[i]) \
                .subs([ (x, midpoints_K_ref[i][0])
                      , (y, midpoints_K_ref[i][1])
                      , (z, midpoints_K_ref[i][2])
                      ]
                     )

  fem_basis = canonical_basis * (A ** -1)
  return fem_basis

v0 = sp.Matrix([   0,   0,   0])
v1 = sp.Matrix([ 0.5,   0,   0])
v2 = sp.Matrix([   1,   0,   0])
v3 = sp.Matrix([   0, 0.5,   0])
v4 = sp.Matrix([ 0.5, 0.5,   0])
v5 = sp.Matrix([   0,   1,   0])
v6 = sp.Matrix([   0,   0, 0.5])
v7 = sp.Matrix([ 0.5,   0, 0.5])
v8 = sp.Matrix([   0, 0.5, 0.5])
v9 = sp.Matrix([   0,   0,   1])

white_up0  = [v0, v1, v3, v6]
white_up1  = [v1, v2, v4, v7]
white_up2  = [v3, v4, v5, v8]
white_up3  = [v6, v7, v8, v9]
green_up   = [v1, v3, v6, v7]
blue_up    = [v1, v3, v4, v7]
blue_down  = [v3, v6, v7, v8]
green_down = [v3, v4, v7, v8]

bwu0 = basis(white_up0)
bwu1 = basis(white_up1)
bwu2 = basis(white_up2)
bwu3 = basis(white_up3)
bgu  = basis(green_up)
bbu  = basis(blue_up)
bbd  = basis(blue_down)
bgd  = basis(green_down)

iwu0 = ((x,   0, 0.5), (y,     0,       0.5-x), (z,       0,       0.5-y-x))
iwu1 = ((x, 0.5,   1), (y,     0, 0.5-(x-0.5)), (z,       0, 0.5-y-(x-0.5)))
iwu2 = ((x,   0, 0.5), (y,   0.5,         1-x), (z,       0, 0.5-(y-0.5)-x))
iwu3 = ((x,   0, 0.5), (y,     0,       0.5-x), (z,     0.5,         1-y-x))
igu  = ((x,   0, 0.5), (y,     0,       0.5-x), (z, 0.5-y-x,         0.5-y))
ibu  = ((x,   0, 0.5), (y, 0.5-x,         0.5), (z,       0,         0.5-y))
ibd  = ((x,   0, 0.5), (y,     0,       0.5-x), (z,   0.5-y,           0.5))
igd  = ((x,   0, 0.5), (y, 0.5-x,         0.5), (z,   0.5-y,         1-x-y))

V = 1/6/8
assert(math.isclose(V, sp.integrate(sp.integrate(sp.integrate(1, iwu0[2]), iwu0[1]), iwu0[0])))
assert(math.isclose(V, sp.integrate(sp.integrate(sp.integrate(1, iwu1[2]), iwu1[1]), iwu1[0])))
assert(math.isclose(V, sp.integrate(sp.integrate(sp.integrate(1, iwu2[2]), iwu2[1]), iwu2[0])))
assert(math.isclose(V, sp.integrate(sp.integrate(sp.integrate(1, iwu3[2]), iwu3[1]), iwu3[0])))
assert(math.isclose(V, sp.integrate(sp.integrate(sp.integrate(1, igu [2]), igu [1]), igu [0])))
assert(math.isclose(V, sp.integrate(sp.integrate(sp.integrate(1, ibu [2]), ibu [1]), ibu [0])))
assert(math.isclose(V, sp.integrate(sp.integrate(sp.integrate(1, ibd [2]), ibd [1]), ibd [0])))
assert(math.isclose(V, sp.integrate(sp.integrate(sp.integrate(1, igd [2]), igd [1]), igd [0])))

tets = [ ("white up 0", white_up0 , bwu0, iwu0)
       , ("white up 1", white_up1 , bwu1, iwu1)
       , ("white up 2", white_up2 , bwu2, iwu2)
       , ("white up 3", white_up3 , bwu3, iwu3)
       , ("green up"  , green_up  , bgu , igu )
       , ("blue up"   , blue_up   , bbu , ibu )
       , ("blue down" , blue_down , bbd , ibd )
       , ("green down", green_down, bgd , igd )
       ]

for desc, verts, basis, bounds in tets:
  print(desc)
  el_mat = sp.zeros(6, 6)

  for i in range(6):
    for j in range(6):
      phi_i = basis.col(i)
      phi_j = basis.col(j)

      val =  sp.integrate(
             sp.integrate(
             sp.integrate( curl(phi_i).dot(curl(phi_j))
                         , bounds[2])
                         , bounds[1])
                         , bounds[0])
      el_mat[i, j] = val

  sp.pprint(el_mat)
  print('')
