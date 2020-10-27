#!/usr/bin/env python3

import numpy as np
import quadpy

# Tested with quadpy 0.16.2
# The syntax changed in more recent versions of quadpy

def quadrature_rule(scheme_name: str, geometry: str, precision=20, comment='') -> str:
  """Returns a string containing C++-formatted quadrature points and weights."""

  if geometry == 'tet':
    scheme = quadpy.t3.schemes[scheme_name]()
    reference_tet = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    point_type = 'Point3D'
  elif geometry == 'triangle':
    scheme = quadpy.t2.schemes[scheme_name]()
    reference_tet = np.asarray([[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]])
    point_type = 'Point2D'
  else:
    raise Exception('Invalid geometry.')

  vol = quadpy.tn.get_vol(reference_tet)
  tet_points = quadpy.tn.transform(scheme.points, reference_tet.T).T
  tet_weights = vol * scheme.weights

  N_int_points = tet_weights.shape[0]

  weight_string = np.array2string(tet_weights, precision=precision, separator=',', max_line_width=1)[1:-1]

  points_string = ""
  for i in range(N_int_points):
    point_string = np.array2string(tet_points[i,:], precision=precision, separator=',', max_line_width=100)[1:-1]
    points_string += f'{point_type}( {{{point_string}}} )'

    if i != N_int_points - 1:
      points_string += ",\n"

  code = ''
  if comment:
    code += "f/// {comment}\n"
  code += f"/// Generated using quadpy, version {quadpy.__version__}\n"
  code += f"/// - scheme:    {scheme.name}\n"
  code += f"/// - key:       {scheme_name}\n"
  code += f"/// - degree:    {scheme.degree}\n"
  code += f"/// - precision: {precision}\n\n"
  code += f'static const std::array< real_t, {N_int_points} > {scheme_name}_weights = {{\n{weight_string}\n}};\n\n'
  code += f'static const std::array< {point_type}, {N_int_points} > {scheme_name}_points = {{\n{points_string}\n}};\n'

  return code


# print(quadrature_rule('hammer_marlowe_stroud_5', 'triangle'))
print(quadrature_rule('lyness_jespersen_15', 'triangle'))
