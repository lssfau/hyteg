#!/usr/bin/env python3

import numpy as np
import quadpy

# Tested with quadpy 0.14.11
# The syntax changed in more recent versions of quadpy

scheme = quadpy.tetrahedron.stroud_t3_5_1()

reference_tet = np.asarray([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])

vol = quadpy.nsimplex.get_vol(reference_tet)
tet_points = quadpy.nsimplex.transform(scheme.points.T, reference_tet.T).T
tet_weights = vol * scheme.weights

N_int_points = tet_weights.shape[0]

weight_string = np.array2string(tet_weights, precision=20, separator=',', max_line_width=1)[1:-1]

points_string = ""
for i in range(N_int_points):
  point_string = np.array2string(tet_points[i,:], precision=20, separator=',', max_line_width=100)[1:-1]
  points_string += "Point3D( {{{}}} )".format(point_string)

  if i != N_int_points - 1:
    points_string += ",\n"

code = "/// 3DT-4: exact for polynomial integrands up to order 5\n"
code += "/// Generated using quadpy 0.14.11 and the tetrahedron.stroud_t3_5_1() scheme\n"
code += "static const std::array< real_t, {} > T4_weights = {{{}}};\n\n".format(N_int_points, weight_string)
code += "static const std::array< Point3D, {} > T4_points = {{{}}};\n".format(N_int_points, points_string)

print(code)