// // #include <array>
// // #include <cmath>
// // #include <iostream>
// // #include <vector>

// // struct Point3D {
// //     double x, y, z;
// // };

// // // Shape-functions N_i(r,s,t) for 13-node quadratic pyramid
// // // and their partial derivatives dN/dr, dN/ds, dN/dt
// // inline void shapeFunctionsQuadratic(double r, double s, double t, double
// // N[13],
// //                                     double dNdr[13], double dNds[13],
// //                                     double dNdt[13]) {
// //     double L = 1 - r - s;

// //     N[0] = L * (1 - t) * (1 - 2 * (r + s + t));
// //     N[1] = r * (1 - t) * (2 * r - 1);
// //     N[2] = s * (1 - t) * (2 * s - 1);
// //     N[3] = t * (2 * t - 1);

// //     N[4] = 4 * r * L * (1 - t);
// //     N[5] = 4 * r * s * (1 - t);
// //     N[6] = 4 * s * L * (1 - t);

// //     N[7] = 4 * r * (1 - t) * t;
// //     N[8] = 4 * L * (1 - t) * t;
// //     N[9] = 4 * s * (1 - t) * t;
// //     N[10] = 4 * r * s * t;
// //     N[11] = 4 * s * L * t;
// //     N[12] = 4 * r * L * t;

// //     // derivatives (computed symbolically)
// //     // dN0/dr
// //     dNdr[0] = -(1 - t) * (1 - 2 * (r + s + t)) - 2 * L * (1 - t);
// //     // dN1/dr
// //     dNdr[1] = (1 - t) * (2 * r - 1) + 2 * r * (1 - t);
// //     // dN2/dr
// //     dNdr[2] = 0;
// //     // dN3/dr
// //     dNdr[3] = 0;
// //     // dN4/dr = 4*(L*(1 - t) - r*(1 - t))
// //     dNdr[4] = 4 * ((-1) * (1 - t) * r + L * (1 - t));
// //     // dN5/dr = 4*(s*(1 - t))
// //     dNdr[5] = 4 * s * (1 - t);
// //     // dN6/dr = 4*0
// //     dNdr[6] = 0;
// //     // dN7/dr = 4*(1 - t)*t
// //     dNdr[7] = 4 * (1 - t) * t;
// //     // dN8/dr = -4*(1 - t)*t
// //     dNdr[8] = -4 * (1 - t) * t;
// //     // dN9/dr = 0
// //     dNdr[9] = 0;
// //     // dN10/dr = 4*s*t
// //     dNdr[10] = 4 * s * t;
// //     // dN11/dr = -4*s*t
// //     dNdr[11] = -4 * s * t;
// // // dN12/dr = 4*L*t - 4*r*t
// // dN12:
// //     dNdr[12] = 4 * ((-1) * t * r + L * t);

// //     // dN/ds similarly (omitted for brevity; compute accordingly)
// //     // For production code, fill all dNds and dNdt entries correctly.
// // }

// // // Invert mapping for quadratic pyramid via Newton iterations
// // bool invertMappingQuadraticPyramid(const std::array<Point3D, 13>& nodes,
// //                                    const Point3D& p, double& r, double& s,
// //                                    double& t, double tol = 1e-8,
// //                                    int maxIter = 20) {
// //     // initial guess: centroid of master
// //     r = s = t = 1.0 / 3.0;
// //     double N[13], dNdr[13], dNds[13], dNdt[13];
// //     for (int iter = 0; iter < maxIter; ++iter) {
// //         shapeFunctionsQuadratic(r, s, t, N, dNdr, dNds, dNdt);
// //         // compute F = X(r,s,t) - p
// //         double Fx = -p.x, Fy = -p.y, Fz = -p.z;
// //         double J[3][3] = {{0}};
// //         for (int i = 0; i < 13; ++i) {
// //             Fx += N[i] * nodes[i].x;
// //             Fy += N[i] * nodes[i].y;
// //             Fz += N[i] * nodes[i].z;
// //             J[0][0] += dNdr[i] * nodes[i].x;
// //             J[0][1] += dNds[i] * nodes[i].x;
// //             J[0][2] += dNdt[i] * nodes[i].x;
// //             J[1][0] += dNdr[i] * nodes[i].y;
// //             J[1][1] += dNds[i] * nodes[i].y;
// //             J[1][2] += dNdt[i] * nodes[i].y;
// //             J[2][0] += dNdr[i] * nodes[i].z;
// //             J[2][1] += dNds[i] * nodes[i].z;
// //             J[2][2] += dNdt[i] * nodes[i].z;
// //         }
// //         double normF = std::sqrt(Fx * Fx + Fy * Fy + Fz * Fz);
// //         if (normF < tol) return true;
// //         // solve J*delta = F
// //         // compute inverse or use Cramer
// //         double detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
// //                       J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
// //                       J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
// //         if (fabs(detJ) < 1e-12) break;
// //         // compute deltas by Cramer's rule
// //         auto detReplace = [&](int col, double Fx, double Fy, double Fz) {
// //             double M[3][3];
// //             for (int i = 0; i < 3; i++)
// //                 for (int j = 0; j < 3; j++)
// //                     M[i][j] = (j == col ? (i == 0   ? Fx
// //                                            : i == 1 ? Fy
// //                                                     : Fz)
// //                                         : J[i][j]);
// //             return M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) -
// //                    M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) +
// //                    M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
// //         };
// //         double dr = detReplace(0, Fx, Fy, Fz) / detJ;
// //         double ds = detReplace(1, Fx, Fy, Fz) / detJ;
// //         double dt = detReplace(2, Fx, Fy, Fz) / detJ;
// //         r -= dr;
// //         s -= ds;
// //         t -= dt;
// //     }
// //     return false;
// // }

// // int main() {
// //     // Example usage: provide 13 nodes in order N0..N12 on master
// //     std::array<Point3D, 13> quadNodes = {{
// //         // physical coordinates for each master node
// //         // fill with your element's 13 points
// //     }};
// //     Point3D p{0.6, 0.4, 0.2};
// //     double r, s, t;
// //     if (invertMappingQuadraticPyramid(quadNodes, p, r, s, t)) {
// //         std::cout << "(r,s,t) = (" << r << "," << s << "," << t << ")\n";
// //     } else {
// //         std::cout << "No convergence" << std::endl;
// //     }
// //     return 0;
// // }
// // // #pragma unroll 4
// // //     for (int i = 0; i < number_of_vertices_of_pyramid; i++) {
// // //         right_part_local[i] = Integrator::integrateForRightPartVector(
// // //             nodes, finite_elements,
// // //             std::make_tuple(index_of_finite_element, i),
// // //             basis_functions_type, basis_functions_element_type);
// // //     }

// // // Point const center_of_base = {
// // //     (node_0.x + node_1.x + node_2.x + node_3.x) / 4.0,
// // //     (node_0.y + node_1.y + node_2.y + node_3.y) / 4.0,
// // //     (node_0.z + node_1.z + node_2.z + node_3.z) / 4.0};

// // // Point const vertice_0_minus_1 = {node_1.x - node_0.x, node_1.y -
// // // node_0.y,
// // //                                  node_1.z - node_0.z};

// // // Point const vertice_0_minus_4 = {node_2.x - node_0.x, node_2.y -
// // // node_0.y,
// // //                                  node_2.z - node_0.z};

// // // Point const vertice_0_minus_point = {point.x - node_0.x, point.y -
// // // node_0.y,
// // //                                      point.z - node_0.z};

// // // double const norm_of_vertice_0_minus_1 =
// // //     vertice_0_minus_1.x * vertice_0_minus_1.x +
// // //     vertice_0_minus_1.y * vertice_0_minus_1.y +
// // //     vertice_0_minus_1.z * vertice_0_minus_1.z;

// // // double const norm_of_vertice_0_minus_4 =
// // //     vertice_0_minus_4.x * vertice_0_minus_4.x +
// // //     vertice_0_minus_4.y * vertice_0_minus_4.y +
// // //     vertice_0_minus_4.z * vertice_0_minus_4.z;

// // // double const ksi_parameter =
// // //     (vertice_0_minus_point.x * vertice_0_minus_1.x +
// // //      vertice_0_minus_point.y * vertice_0_minus_1.y +
// // //      vertice_0_minus_point.z * vertice_0_minus_1.z) /
// // //     norm_of_vertice_0_minus_1;

// // // double const nu_parameter =
// // //     (vertice_0_minus_point.x * vertice_0_minus_4.x +
// // //      vertice_0_minus_point.y * vertice_0_minus_4.y +
// // //      vertice_0_minus_point.z * vertice_0_minus_4.z) /
// // //     norm_of_vertice_0_minus_4;

// // // Point const vetice_vertex_minus_center = {(vertex.x -
// center_of_base.x),
// // //                                           (vertex.y -
// center_of_base.y),
// // //                                           (vertex.z -
// center_of_base.z)};

// // // double const norm_of_vetice_vertex_minus_center =
// // //     vetice_vertex_minus_center.x * vetice_vertex_minus_center.x +
// // //     vetice_vertex_minus_center.y * vetice_vertex_minus_center.y +
// // //     vetice_vertex_minus_center.z * vetice_vertex_minus_center.z;

// // // Point const vetice_point_minus_center = {point.x - center_of_base.x,
// // //                                          point.y - center_of_base.y,
// // //                                          point.z - center_of_base.z};

// // // double const tetta_parameter =
// // //     (vetice_point_minus_center.x * vetice_vertex_minus_center.x +
// // //      vetice_point_minus_center.y * vetice_vertex_minus_center.y +
// // //      vetice_point_minus_center.z * vetice_vertex_minus_center.z) /
// // //     norm_of_vetice_vertex_minus_center;

// // void shapeFunctionsQuadratic1(double r, double s, double t, double N[13],
// //                               double dNdr[13], double dNds[13],
// //                               double dNdt[13]) {
// //     // N[0] = (1 - r) * (1 - s) * (1 - t);
// //     // N[1] = r * (1 - s) * (1 - t);
// //     // N[2] = (1 - r) * s * (1 - t);
// //     // N[3] = r * s * (1 - t);
// //     // N[4] = t;

// //     // dNdr[0] = -(1 - s) * (1 - t);
// //     // dNdr[1] = (1 - s) * (1 - t);
// //     // dNdr[2] = -s * (1 - t);
// //     // dNdr[3] = s * (1 - t);
// //     // dNdr[4] = 0;

// //     // dNds[0] = -(1 - r) * (1 - t);
// //     // dNds[1] = -r * (1 - t);
// //     // dNds[2] = (1 - r) * (1 - t);
// //     // dNds[3] = r * (1 - t);
// //     // dNds[4] = 0;

// //     // dNdt[0] = -(1 - s) * (1 - r);
// //     // dNdt[1] = -(1 - s) * r;
// //     // dNdt[2] = -(1 - r) * s;
// //     // dNdt[3] = -r * s;
// //     // dNdt[4] = 1;

// //     N[0] = 8 * (0.5 - r) * (1 - r) * (0.5 - s) * (1 - s) * (0.5 - t) * (1
// -
// //     t); N[1] = 16 * r * (1 - r) * (0.5 - s) * (1 - s) * (0.5 - t) * (1 -
// t);
// //     N[2] = -8 * r * (0.5 - r) * (0.5 - s) * (1 - s) * (0.5 - t) * (1 - t);
// //     N[3] = 16 * (0.5 - r) * (1 - r) * s * (1 - s) * (0.5 - t) * (1 - t);
// //     N[4] = -16 * r * (0.5 - r) * s * (1 - s) * (0.5 - t) * (1 - t);
// //     N[5] = -8 * (0.5 - r) * (1 - r) * s * (0.5 - s) * (0.5 - t) * (1 - t);
// //     N[6] = -16 * r * (1 - r) * s * (0.5 - s) * (0.5 - t) * (1 - t);
// //     N[7] = 8 * r * (0.5 - r) * s * (0.5 - s) * (0.5 - t) * (1 - t);
// //     N[8] = 16 * t * (1 - t) * (0.75 - r) * (0.75 - s);
// //     N[9] = -16 * t * (1 - t) * (0.25 - r) * (0.75 - s);
// //     N[10] = -16 * t * (1 - t) * (0.75 - r) * (0.25 - s);
// //     N[11] = 16 * t * (1 - t) * (0.25 - r) * (0.25 - s);
// //     N[12] = -2 * t * (0.5 - t);

// //     dNdr[0] = (-12 + 16 * r) * (0.5 - s) * (1 - s) * (0.5 - t) * (1 - t);
// //     dNdr[1] = (16 - 32 * r) * (0.5 - s) * (1 - s) * (0.5 - t) * (1 - t);
// //     dNdr[2] = (-4 + 16 * r) * (0.5 - s) * (1 - s) * (0.5 - t) * (1 - t);
// //     dNdr[3] = (-24 + 32 * r) * s * (1 - s) * (0.5 - t) * (1 - t);
// //     dNdr[4] = (-8 + 32 * r) * s * (1 - s) * (0.5 - t) * (1 - t);
// //     dNdr[5] = (12 - 16 * r) * s * (0.5 - s) * (0.5 - t) * (1 - t);
// //     dNdr[6] = (-16 + 32 * r) * s * (0.5 - s) * (0.5 - t) * (1 - t);
// //     dNdr[7] = (4 - 16 * r) * s * (0.5 - s) * (0.5 - t) * (1 - t);
// //     dNdr[8] = -16 * t * (1 - t) * (0.75 - s);
// //     dNdr[9] = 16 * t * (1 - t) * (0.75 - s);
// //     dNdr[10] = 16 * t * (1 - t) * (0.25 - s);
// //     dNdr[11] = -16 * t * (1 - t) * (0.25 - s);
// //     dNdr[12] = 0;

// //     dNds[0] = (0.5 - r) * (1 - r) * (-12 + 16 * s) * (0.5 - t) * (1 - t);
// //     dNds[1] = r * (1 - r) * (-24 + 32 * s) * (0.5 - t) * (1 - t);
// //     dNds[2] = r * (0.5 - r) * (12 - 16 * s) * (0.5 - t) * (1 - t);
// //     dNds[3] = (0.5 - r) * (1 - r) * (16 - 32 * s) * (0.5 - t) * (1 - t);
// //     dNds[4] = r * (0.5 - r) * (-16 + 32 * s) * (0.5 - t) * (1 - t);
// //     dNds[5] = (0.5 - r) * (1 - r) * (-4 + 16 * s) * (0.5 - t) * (1 - t);
// //     dNds[6] = r * (1 - r) * (-8 + 32 * s) * (0.5 - t) * (1 - t);
// //     dNds[7] = r * (0.5 - r) * (4 - 16 * s) * (0.5 - t) * (1 - t);
// //     dNds[8] = -16 * t * (1 - t) * (0.75 - r);
// //     dNds[9] = 16 * t * (1 - t) * (0.25 - r);
// //     dNds[10] = 16 * t * (1 - t) * (0.75 - r);
// //     dNds[11] = -16 * t * (1 - t) * (0.25 - r);
// //     dNds[12] = 0;

// //     dNdt[0] = (0.5 - r) * (1 - r) * (0.5 - s) * (1 - s) * (-12 + 16 * t);
// //     dNdt[1] = r * (1 - r) * (0.5 - s) * (1 - s) * (-24 + 32 * t);
// //     dNdt[2] = r * (0.5 - r) * (0.5 - s) * (1 - s) * (12 - 16 * t);
// //     dNdt[3] = (0.5 - r) * (1 - r) * s * (1 - s) * (-24 + 32 * t);
// //     dNdt[4] = r * (0.5 - r) * s * (1 - s) * (24 - 32 * t);
// //     dNdt[5] = (0.5 - r) * (1 - r) * s * (0.5 - s) * (12 - 16 * t);
// //     dNdt[6] = r * (1 - r) * s * (0.5 - s) * (24 - 32 * t);
// //     dNdt[7] = r * (0.5 - r) * s * (0.5 - s) * (-12 + 16 * t);
// //     dNdt[8] = (16 - 32 * t) * (0.75 - r) * (0.75 - s);
// //     dNdt[9] = (-16 + 32 * t) * (0.25 - r) * (0.75 - s);
// //     dNdt[10] = (-16 + 32 * t) * (0.75 - r) * (0.25 - s);
// //     dNdt[11] = (16 - 32 * t) * (0.25 - r) * (0.25 - s);
// //     dNdt[12] = -1 + 4 * t;
// // }

// // void shapeFunctionsQuadratic1(double r, double s, double t, double N[13],
// //                               double dNdr[13], double dNds[13],
// //                               double dNdt[13]) {
// //     double R0 = (1 - 2 * r) * (1 - r);
// //     double R1 = 4 * r * (1 - r);
// //     double R2 = -r * (1 - 2 * r);
// //     double dR0 = -3 + 4 * r;
// //     double dR1 = 4 - 8 * r;
// //     double dR2 = -1 + 4 * r;

// //     double S0 = (1 - 2 * s) * (1 - s);
// //     double S1 = 4 * s * (1 - s);
// //     double S2 = -s * (1 - 2 * s);
// //     double dS0 = -3 + 4 * s;
// //     double dS1 = 4 - 8 * s;
// //     double dS2 = -1 + 4 * s;

// //     double T0 = (1 - 2 * t) * (1 - t);
// //     double T1 = 4 * t * (1 - t);
// //     double T2 = -t * (1 - 2 * t);
// //     double dT0 = -3 + 4 * t;
// //     double dT1 = 4 - 8 * t;
// //     double dT2 = -1 + 4 * t;

// //     N[0] = R0 * S0 * T0;
// //     N[1] = R1 * S0 * T0;
// //     N[2] = R2 * S0 * T0;
// //     N[3] = R0 * S1 * T0;
// //     N[4] = R2 * S1 * T0;
// //     N[5] = R0 * S2 * T0;
// //     N[6] = R1 * S2 * T0;
// //     N[7] = R2 * S2 * T0;
// //     N[8] = R0 * S0 * T1;
// //     N[9] = R2 * S0 * T1;
// //     N[10] = R0 * S2 * T1;
// //     N[11] = R2 * S2 * T1;
// //     N[12] = T2;

// //     dNdr[0] = dR0 * S0 * T0;
// //     dNdr[1] = dR1 * S0 * T0;
// //     dNdr[2] = dR2 * S0 * T0;
// //     dNdr[3] = dR0 * S1 * T0;
// //     dNdr[4] = dR2 * S1 * T0;
// //     dNdr[5] = dR0 * S2 * T0;
// //     dNdr[6] = dR1 * S2 * T0;
// //     dNdr[7] = dR2 * S2 * T0;
// //     dNdr[8] = dR0 * S0 * T1;
// //     dNdr[9] = dR2 * S0 * T1;
// //     dNdr[10] = dR0 * S2 * T1;
// //     dNdr[11] = dR2 * S2 * T1;
// //     dNdr[12] = 0;

// //     dNds[0] = R0 * dS0 * T0;
// //     dNds[1] = R1 * dS0 * T0;
// //     dNds[2] = R2 * dS0 * T0;
// //     dNds[3] = R0 * dS1 * T0;
// //     dNds[4] = R2 * dS1 * T0;
// //     dNds[5] = R0 * dS2 * T0;
// //     dNds[6] = R1 * dS2 * T0;
// //     dNds[7] = R2 * dS2 * T0;
// //     dNds[8] = R0 * dS0 * T1;
// //     dNds[9] = R2 * dS0 * T1;
// //     dNds[10] = R0 * dS2 * T1;
// //     dNds[11] = R2 * dS2 * T1;
// //     dNds[12] = 0;

// //     dNdt[0] = R0 * S0 * dT0;
// //     dNdt[1] = R1 * S0 * dT0;
// //     dNdt[2] = R2 * S0 * dT0;
// //     dNdt[3] = R0 * S1 * dT0;
// //     dNdt[4] = R2 * S1 * dT0;
// //     dNdt[5] = R0 * S2 * dT0;
// //     dNdt[6] = R1 * S2 * dT0;
// //     dNdt[7] = R2 * S2 * dT0;
// //     dNdt[8] = R0 * S0 * dT1;
// //     dNdt[9] = R2 * S0 * dT1;
// //     dNdt[10] = R0 * S2 * dT1;
// //     dNdt[11] = R2 * S2 * dT1;
// //     dNdt[12] = dT2;
// // }

// // switch (number_basis_function) {
// //     case 0:
// //         return 8 * (0.5 - ksi_parameter) * (1 - ksi_parameter) *
// //                (0.5 - nu_parameter) * (1 - nu_parameter) *
// //                (0.5 - tetta_parameter) * (1 - tetta_parameter);
// //         break;

// //     case 1:
// //         return 16 * ksi_parameter * (1 - ksi_parameter) * (0.5 -
// //         nu_parameter) *
// //                (1 - nu_parameter) * (0.5 - tetta_parameter) *
// //                (1 - tetta_parameter);
// //         break;

// //     case 2:
// //         return -8 * ksi_parameter * (0.5 - ksi_parameter) *
// //                (0.5 - nu_parameter) * (1 - nu_parameter) *
// //                (0.5 - tetta_parameter) * (1 - tetta_parameter);
// //         break;

// //     case 3:
// //         return 16 * (0.5 - ksi_parameter) * (1 - ksi_parameter) *
// //         nu_parameter *
// //                (1 - nu_parameter) * (0.5 - tetta_parameter) *
// //                (1 - tetta_parameter);
// //         break;

// //     case 4:
// //         return -16 * ksi_parameter * (0.5 - ksi_parameter) * nu_parameter
// *
// //                (1 - nu_parameter) * (0.5 - tetta_parameter) *
// //                (1 - tetta_parameter);
// //         break;

// //     case 5:
// //         return -8 * (0.5 - ksi_parameter) * (1 - ksi_parameter) *
// //         nu_parameter *
// //                (0.5 - nu_parameter) * (0.5 - tetta_parameter) *
// //                (1 - tetta_parameter);
// //         break;

// //     case 6:
// //         return -16 * ksi_parameter * (1 - ksi_parameter) * nu_parameter *
// //                (0.5 - nu_parameter) * (0.5 - tetta_parameter) *
// //                (1 - tetta_parameter);
// //         break;

// //     case 7:
// //         return 8 * ksi_parameter * (0.5 - ksi_parameter) * nu_parameter *
// //                (0.5 - nu_parameter) * (0.5 - tetta_parameter) *
// //                (1 - tetta_parameter);
// //         break;

// //     case 8:
// //         return 16 * tetta_parameter * (1 - tetta_parameter) *
// //                (0.75 - ksi_parameter) * (0.75 - nu_parameter);
// //         break;

// //     case 9:
// //         return -16 * tetta_parameter * (1 - tetta_parameter) *
// //                (0.25 - ksi_parameter) * (0.75 - nu_parameter);
// //         break;

// //     case 10:
// //         return -16 * tetta_parameter * (1 - tetta_parameter) *
// //                (0.75 - ksi_parameter) * (0.25 - nu_parameter);
// //         break;

// //     case 11:
// //         return 16 * tetta_parameter * (1 - tetta_parameter) *
// //                (0.25 - ksi_parameter) * (0.25 - nu_parameter);
// //         break;

// //     case 12:
// //         return -2 * tetta_parameter * (0.5 - tetta_parameter);
// //         break;

// //     default:
// //         return 0.0;
// // }

// // std::vector<double> solveSystem3x3(std::vector<std::vector<double>>& J,
// //                                    std::vector<double>& b) {
// //     std::vector<std::vector<double>> A = J;
// //     std::vector<double> x(3, 0.0);

// //     for (int col = 0; col < 3; ++col) {
// //         int max_row = col;
// //         for (int row = col + 1; row < 3; ++row) {
// //             if (abs(A[row][col]) > abs(A[max_row][col])) {
// //                 max_row = row;
// //             }
// //         }

// //         std::swap(A[col], A[max_row]);
// //         std::swap(b[col], b[max_row]);

// //         double pivot = A[col][col];

// //         if (abs(pivot) < 1e-10) {
// //             std::cout
// //                 << "Матрица вырождена (нет решения или их бесконечно
// много)!"
// //                 << '\n';
// //             return x;
// //         }

// //         for (int k = col; k < 3; ++k) {
// //             A[col][k] /= pivot;
// //         }

// //         b[col] /= pivot;

// //         for (int row = 0; row < 3; ++row) {
// //             if (row != col && abs(A[row][col]) > 1e-10) {
// //                 double factor = A[row][col];

// //                 for (int k = col; k < 3; ++k) {
// //                     A[row][k] -= factor * A[col][k];
// //                 }

// //                 b[row] -= factor * b[col];
// //             }
// //         }
// //     }

// //     x = b;
// //     return x;
// // }
// // bool invertMappingQuadraticPyramid1(const std::array<Point, 13>& nodes,
// //                                     const Point& p, double& r, double& s,
// //                                     double& t, double tol = 1e-8,
// //                                     int maxIter = 20) {
// //     r = s = t = 1.0 / 3.0;

// //     double N[13];
// //     double dNdr[13];
// //     double dNds[13];
// //     double dNdt[13];

// //     for (int iter = 0; iter < 40; ++iter) {
// //         shapeFunctionsQuadratic1(r, s, t, N, dNdr, dNds, dNdt);

// //         double Fx = 0.0;
// //         double Fy = 0.0;
// //         double Fz = 0.0;

// //         std::vector<std::vector<double>> J = {{0, 0, 0}, {0, 0, 0}, {0, 0,
// //         0}};

// //         for (int i = 0; i < 13; ++i) {
// //             Fx += N[i] * nodes[i].x;
// //             Fy += N[i] * nodes[i].y;
// //             Fz += N[i] * nodes[i].z;

// //             J[0][0] += dNdr[i] * nodes[i].x;
// //             J[0][1] += dNds[i] * nodes[i].x;
// //             J[0][2] += dNdt[i] * nodes[i].x;

// //             J[1][0] += dNdr[i] * nodes[i].y;
// //             J[1][1] += dNds[i] * nodes[i].y;
// //             J[1][2] += dNdt[i] * nodes[i].y;

// //             J[2][0] += dNdr[i] * nodes[i].z;
// //             J[2][1] += dNds[i] * nodes[i].z;
// //             J[2][2] += dNdt[i] * nodes[i].z;
// //         }

// //         std::vector<double> d = {p.x - Fx, p.y - Fy, p.z - Fz};

// //         if (std::sqrt(std::pow(d[0], 2) + std::pow(d[1], 2) +
// //                       std::pow(d[2], 2)) < tol) {
// //             return true;
// //         }

// //         std::vector<double> delta = solveSystem3x3(J, d);

// //         r += delta[0];
// //         s += delta[1];
// //         t += delta[2];
// //     }
// //     return false;
// // }

// // // N[0] = (1 - r) * (1 - s) * (1 - t) * (1 - 2 * r - 2 * s);
// // // N[1] = 4 * r * (1 - r) * (1 - s) * (1 - t);
// // // N[2] = r * (1 - s) * (1 - t) * (2 * r - 1 - 2 * s);
// // // N[3] = 4 * (1 - r) * s * (1 - s) * (1 - t);
// // // N[4] = 4 * r * s * (1 - s) * (1 - t);
// // // N[5] = (1 - r) * s * (1 - t) * (2 * s - 1 - 2 * r);
// // // N[6] = 4 * r * (1 - r) * s * (1 - t);
// // // N[7] = r * s * (1 - t) * (2 * r + 2 * s - 3);
// // // N[8] = 4 * (1 - r) * (1 - s) * t * (1 - t);
// // // N[9] = 4 * r * (1 - s) * t * (1 - t);
// // // N[10] = 4 * (1 - r) * s * t * (1 - t);
// // // N[11] = 4 * r * s * t * (1 - t);
// // // N[12] = t * (2 * t - 1);

// // // dNdr[0] = (1 - s) * (1 - t) * (-3 + 4 * r + 2 * s);
// // // dNdr[1] = 4 * (1 - 2 * r) * (1 - s) * (1 - t);
// // // dNdr[2] = (1 - s) * (1 - t) * (4 * r - 1 - 2 * s);
// // // dNdr[3] = -4 * s * (1 - s) * (1 - t);
// // // dNdr[4] = 4 * s * (1 - s) * (1 - t);
// // // dNdr[5] = s * (1 - t) * (-2 * s - 1 + 4 * r);
// // // dNdr[6] = 4 * (1 - 2 * r) * s * (1 - t);
// // // dNdr[7] = s * (1 - t) * (4 * r + 2 * s - 3);
// // // dNdr[8] = -4 * (1 - s) * t * (1 - t);
// // // dNdr[9] = 4 * (1 - s) * t * (1 - t);
// // // dNdr[10] = -4 * s * t * (1 - t);
// // // dNdr[11] = 4 * s * t * (1 - t);
// // // dNdr[12] = 0;

// // // dNds[0] = (1 - r) * (1 - t) * (-3 + 2 * r + 4 * s);
// // // dNds[1] = -4 * r * (1 - r) * (1 - t);
// // // dNds[2] = r * (1 - t) * (-1 - 2 * r + 4 * s);
// // // dNds[3] = 4 * (1 - r) * (1 - 2 * s) * (1 - t);
// // // dNds[4] = 4 * r * (1 - 2 * s) * (1 - t);
// // // dNds[5] = (1 - r) * (1 - t) * (4 * s - 1 - 2 * r);
// // // dNds[6] = 4 * r * (1 - r) * (1 - t);
// // // dNds[7] = r * (1 - t) * (2 * r + 4 * s - 3);
// // // dNds[8] = -4 * (1 - r) * t * (1 - t);
// // // dNds[9] = -4 * r * t * (1 - t);
// // // dNds[10] = 4 * (1 - r) * t * (1 - t);
// // // dNds[11] = 4 * r * t * (1 - t);
// // // dNds[12] = 0;

// // // dNdt[0] = -(1 - r) * (1 - s) * (1 - 2 * r - 2 * s);
// // // dNdt[1] = -4 * r * (1 - r) * (1 - s);
// // // dNdt[2] = -r * (1 - s) * (2 * r - 1 - 2 * s);
// // // dNdt[3] = -4 * (1 - r) * s * (1 - s);
// // // dNdt[4] = -4 * r * s * (1 - s);
// // // dNdt[5] = -(1 - r) * s * (2 * s - 1 - 2 * r);
// // // dNdt[6] = -4 * r * (1 - r) * s;
// // // dNdt[7] = -r * s * (2 * r + 2 * s - 3);
// // // dNdt[8] = 4 * (1 - r) * (1 - s) * (1 - 2 * t);
// // // dNdt[9] = 4 * r * (1 - s) * (1 - 2 * t);
// // // dNdt[10] = 4 * (1 - r) * s * (1 - 2 * t);
// // // dNdt[11] = 4 * r * s * (1 - 2 * t);
// // // dNdt[12] = 4 * t - 1;

// // double right_part_test(double ksi, double nu, double tetta,
// //                        std::array<Point, 5> pyramide_nodes) {
// //     double x = 0.0;
// //     double y = 0.0;
// //     double z = 0.0;

// //     Point nodec = {(pyramide_nodes[0].x + pyramide_nodes[1].x +
// //                     pyramide_nodes[2].x + pyramide_nodes[3].x) /
// //                        4.0,
// //                    (pyramide_nodes[0].y + pyramide_nodes[1].y +
// //                     pyramide_nodes[2].y + pyramide_nodes[3].y) /
// //                        4.0,
// //                    (pyramide_nodes[0].z + pyramide_nodes[1].z +
// //                     pyramide_nodes[2].z + pyramide_nodes[3].z) /
// //                        4.0};

// //     double D10 = pow(pyramide_nodes[1].x - pyramide_nodes[0].x, 2) +
// //                  pow(pyramide_nodes[1].y - pyramide_nodes[0].y, 2) +
// //                  pow(pyramide_nodes[1].z - pyramide_nodes[0].z, 2);

// //     double D20 = pow(pyramide_nodes[2].x - pyramide_nodes[0].x, 2) +
// //                  pow(pyramide_nodes[2].y - pyramide_nodes[0].y, 2) +
// //                  pow(pyramide_nodes[2].z - pyramide_nodes[0].z, 2);

// //     double D4c = pow(pyramide_nodes[4].x - nodec.x, 2) +
// //                  pow(pyramide_nodes[4].y - nodec.y, 2) +
// //                  pow(pyramide_nodes[4].z - nodec.z, 2);

// //     std::vector<std::vector<double>> A(3, std::vector<double>(3));

// //     A[0][0] = (pyramide_nodes[1].x - pyramide_nodes[0].x) / D10;
// //     A[0][1] = (pyramide_nodes[1].y - pyramide_nodes[0].y) / D10;
// //     A[0][2] = (pyramide_nodes[1].z - pyramide_nodes[0].z) / D10;

// //     A[1][0] = (pyramide_nodes[2].x - pyramide_nodes[0].x) / D20;
// //     A[1][1] = (pyramide_nodes[2].y - pyramide_nodes[0].y) / D20;
// //     A[1][2] = (pyramide_nodes[2].z - pyramide_nodes[0].z) / D20;

// //     A[2][0] = (pyramide_nodes[4].x - nodec.x) / D4c;
// //     A[2][1] = (pyramide_nodes[4].y - nodec.y) / D4c;
// //     A[2][2] = (pyramide_nodes[4].z - nodec.z) / D4c;

// //     std::vector<double> b(3);

// //     b[0] = ksi + pyramide_nodes[0].x * A[0][0] + pyramide_nodes[0].y *
// //     A[0][1] +
// //            pyramide_nodes[0].z * A[0][2];
// //     b[1] = nu + pyramide_nodes[0].x * A[1][0] + pyramide_nodes[0].y *
// A[1][1]
// //     +
// //            pyramide_nodes[0].z * A[1][2];
// //     b[2] = tetta + nodec.x * A[2][0] + nodec.y * A[2][1] + nodec.z *
// A[2][2];

// //     double first_determinant_component = A[0][0] * A[1][1] * A[2][2];
// //     double second_determinant_component = A[0][1] * A[1][2] * A[2][0];
// //     double third_determinant_component = A[0][2] * A[1][1] * A[2][0];
// //     double fourth_determinant_component = A[0][0] * A[1][2] * A[2][1];
// //     double fifth_determinant_component = A[0][2] * A[2][1] * A[1][0];
// //     double sixth_determinant_component = A[0][1] * A[2][2] * A[1][0];

// //     double determinant =
// //         (first_determinant_component + second_determinant_component -
// //          third_determinant_component - fourth_determinant_component +
// //          fifth_determinant_component - sixth_determinant_component);

// //     std::vector<std::vector<double>> A_Cofactor(3,
// std::vector<double>(3));

// //     A_Cofactor[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]);
// //     A_Cofactor[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]);
// //     A_Cofactor[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]);

// //     A_Cofactor[1][0] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]);
// //     A_Cofactor[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]);
// //     A_Cofactor[1][2] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]);

// //     A_Cofactor[2][0] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
// //     A_Cofactor[2][1] = -(A[0][0] * A[2][1] - A[0][1] * A[2][0]);
// //     A_Cofactor[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]);

// //     std::vector<std::vector<double>> A_inv(3, std::vector<double>(3));

// //     for (int i = 0; i < 3; i++) {
// //         for (int j = 0; j < 3; j++) {
// //             A_inv[i][j] = A_Cofactor[j][i] / determinant;
// //         }
// //     }

// //     x = A_inv[0][0] * b[0] + A_inv[0][1] * b[1] + A_inv[0][2] * b[2];
// //     y = A_inv[1][0] * b[0] + A_inv[1][1] * b[1] + A_inv[1][2] * b[2];
// //     z = A_inv[2][0] * b[0] + A_inv[2][1] * b[1] + A_inv[2][2] * b[2];

// //     return MathUtils::calculateF({x, y, z});
// // }

// // double computeJacobianPhysicalToLocal(Point node_0, Point node_1, Point
// // node_2,
// //     Point node_3, Point node_4) {
// // const Point node_c = {(node_0.x + node_1.x + node_2.x + node_3.x) / 4.0,
// // (node_0.y + node_1.y + node_2.y + node_3.y) / 4.0,
// // (node_0.z + node_1.z + node_2.z + node_3.z) / 4.0};

// // const double D10 = pow(node_1.x - node_0.x, 2) +
// // pow(node_1.y - node_0.y, 2) +
// // pow(node_1.z - node_0.z, 2);

// // const double D20 = pow(node_2.x - node_0.x, 2) +
// // pow(node_2.y - node_0.y, 2) +
// // pow(node_2.z - node_0.z, 2);

// // const double D4c = pow(node_4.x - node_c.x, 2) +
// // pow(node_4.y - node_c.y, 2) +
// // pow(node_4.z - node_c.z, 2);

// // const Point minus10 = {node_1.x - node_0.x, node_1.y - node_0.y,
// // node_1.z - node_0.z};

// // const Point minus20 = {node_2.x - node_0.x, node_2.y - node_0.y,
// // node_2.z - node_0.z};

// // const Point minus4c = {node_4.x - node_c.x, node_4.y - node_c.y,
// // node_4.z - node_c.z};

// // std::vector<std::vector<double>> J(3, std::vector<double>(3));

// // J[0][0] = D10 / ((minus10.x == 0) ? 1 : minus10.x);
// // J[0][1] = D20 / ((minus20.x == 0) ? 1 : minus20.x);
// // J[0][2] = D4c / ((minus4c.x == 0) ? 1 : minus4c.x);

// // J[1][0] = D10 / ((minus10.y == 0) ? 1 : minus10.y);
// // J[1][1] = D20 / ((minus20.y == 0) ? 1 : minus20.y);
// // J[1][2] = D4c / ((minus4c.y == 0) ? 1 : minus4c.y);

// // J[2][0] = D10 / ((minus10.z == 0) ? 1 : minus10.z);
// // J[2][1] = D20 / ((minus20.z == 0) ? 1 : minus20.z);
// // J[2][2] = D4c / ((minus4c.z == 0) ? 1 : minus4c.z);

// // return J[0][0] * (J[1][1] * J[2][2] - J[2][1] * J[1][2]) -
// // J[1][0] * (J[0][1] * J[2][2] - J[2][1] * J[0][2]) +
// // J[2][0] * (J[0][1] * J[1][2] - J[1][1] * J[0][2]);
// // }

// #include <algorithm>
// #include <iostream>
// #include <set>
// #include <vector>

// using namespace std;

// // 0 1 2 3
// // 2 3 5 6
// // 5 6 8 9
// // 3 4 6 7
// // 6 7 9 10

// int main() {
//     int slae_size = 11;
//     std::vector<int> ig(11), jg;
//     int lower_triangle_size = 0;
//     vector<set<int>> connections = {
//         {1, 2, 3},     {2, 3},  {3, 5, 6}, {4, 5, 6, 7}, {6, 7}, {6, 8, 9},
//         {7, 8, 9, 10}, {9, 10}, {9},       {10},         {}};

//     for (int i = 0; i < 11; i++) {
//         lower_triangle_size += connections[i].size();
//     }
//     jg.resize(lower_triangle_size);
//     int current_column_index = 0;

//     ig[0] = 0;

//     for (int i = 0; i < slae_size; i++) {
//         int current_count_of_elements = 0;

// #pragma unroll 4
//         for (int j = 0; j <= i; j++) {
//             if (count(connections[j].begin(), connections[j].end(), i) != 0)
//             {
//                 jg[current_column_index] = j;
//                 current_column_index++;
//                 current_count_of_elements++;
//             }
//         }

//         ig[i + 1] = ig[i] + current_count_of_elements;
//     }

//     // cout << jg.size() << '\n';
//     // for (auto index : ig) {
//     //     cout << index << ' ';
//     // }
//     // cout << '\n';
//     // for (auto index : jg) {
//     //     cout << index << ' ';
//     // }
//     // cout << '\n';
//     // for (int i = 0; i < ig.size(); i++) {
//     //     cout << i << " " << ig[i + 1] - ig[i] << '\n';
//     // }
//     return 0;
// }

//     // J[0][0] = (node_1.x - node_0.x) / D10;  // d_ksi_d_x
//     // J[0][1] = (node_2.x - node_0.x) / D20;  // d_nu_d_x
//     // J[0][2] = (node_4.x - node_c.x) / D4c;  // d_tetta_d_x

//     // J[1][0] = (node_1.y - node_0.y) / D10;  // d_ksi_d_y
//     // J[1][1] = (node_2.y - node_0.y) / D20;  // d_nu_d_y
//     // J[1][2] = (node_4.y - node_c.y) / D4c;  // d_tetta_d_y

//     // J[2][0] = (node_1.z - node_0.z) / D10;  // d_ksi_d_z
//     // J[2][1] = (node_2.z - node_0.z) / D20;  // d_nu_d_z
//     // J[2][2] = (node_4.z - node_c.z) / D4c;  // d_tetta_d_z

//     std::vector<Point> nodes_0 = {
//         {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0.5, 0.5, 0.5}};

//     std::vector<Point> nodes_1 = {
//         {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1}, {0.5, 0.5, 0.5}};

//     std::vector<Point> nodes_2 = {
//         {0, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 1}, {0.5, 0.5, 0.5}};

//     std::vector<Point> nodes_3 = {
//         {1, 0, 0}, {1, 1, 0}, {1, 0, 1}, {1, 1, 1}, {0.5, 0.5, 0.5}};

//     std::vector<Point> points = {
//         {0, 0, 0},         {0.5, 0, 0},       {1, 0, 0},
//         {0, 0.5, 0},       {1, 0.5, 0},       {0, 1, 0},
//         {0.5, 1, 0},       {1, 1, 0},         {0.25, 0.25, 0.5},
//         {0.75, 0.25, 0.5}, {0.25, 0.75, 0.5}, {0.75, 0.75, 0.5},
//         {0.5, 0.5, 1}};

//     // std::vector<double> result = {
//     //     0,      0.25,   1,      0.25,   0.5,    1.25,   1,      1.25, 2,
//     //     0.1875, 0.6875, 0.6875, 1.1875, 0.25,   0.5,    1.25,   0.5,  1.5,
//     //     1.25,   1.5,    2.25,   0.6875, 1.1875, 1.1875, 1.6875,
//     1,    1.25,
//     //     2,      1.25,   1.5,    2.25,   2,      2.25,   3,      0.75};
//     std::vector<double> result = {0,      0.25,   1,      0.25, 0.5,
//                                   1.25,   1,      1.25,   2,    0.1875,
//                                   0.6875, 0.6875, 1.1875, 0.75};

//     std::vector<int> indexes_0 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12,
//     13}; std::vector<int> indexes_1 = {4, 5, 6, 7, 8}; std::vector<int>
//     indexes_2 = {0, 2, 4, 6, 8}; std::vector<int> indexes_3 = {1, 3, 5, 7,
//     8};

//     double ksi = (test_points[0].x - 0) / (1 - 0);
//     double nu = (test_points[0].y - 0) / (1 - 0);
//     double tetta = (test_points[0].z - 0) / (0.5 - 0);
//     // std::cout << " " << ksi << " " << nu << " " << tetta << '\n';
//     double ksi_n = ksi / (1 - tetta);
//     double nu_n = nu / (1 - tetta);
//     double tetta_n = tetta;
//     // std::cout << " " << ksi_n << " " << nu_n << " " << tetta_n << '\n';
//     double sum = 0;
//     for (int i = 0; i < 14; i++) {
//         // sum += getLinearBasisFunction11(test_points[0], nodes_0[0],
//         // nodes_0[1],
//         //                                 nodes_0[2], nodes_0[3],
//         nodes_0[4],
//         //                                 i);

//         sum += result[indexes_0[i]] *
//                getLinearBasisFunction11(test_points[0], nodes_0[0],
//                nodes_0[1],
//                                         nodes_0[2], nodes_0[3], nodes_0[4],
//                                         i);
//     }
//     std::cout << sum << " "
//               << test_points[0].x * test_points[0].x +
//                      test_points[0].y * test_points[0].y +
//                      test_points[0].z * test_points[0].z
//               << '\n';
//     // sum = 0;
//     // for (int i = 0; i < 5; i++) {
//     //     sum += result[indexes_1[i]] *
//     //            getLinearBasisFunction11(test_points[1], nodes_1[0],
//     //            nodes_1[1],
//     //                                     nodes_1[2], nodes_1[3],
//     nodes_1[4],
//     //                                     i);
//     // }
//     // std::cout << sum << " "
//     //           << test_points[1].x + test_points[1].y + test_points[1].z <<
//     //           '\n';
//     // sum = 0;
//     // for (int i = 0; i < 5; i++) {
//     //     sum += result[indexes_2[i]] *
//     //            getLinearBasisFunction11(test_points[2], nodes_2[0],
//     //            nodes_2[1],
//     //                                     nodes_2[2], nodes_2[3],
//     nodes_2[4],
//     //                                     i);
//     // }
//     // std::cout << sum << " "
//     //           << test_points[2].x + test_points[2].y + test_points[2].z <<
//     //           '\n';
//     // sum = 0;
//     // for (int i = 0; i < 5; i++) {
//     //     sum += result[indexes_3[i]] *
//     //            getLinearBasisFunction11(test_points[3], nodes_3[0],
//     //            nodes_3[1],
//     //                                     nodes_3[2], nodes_3[3],
//     nodes_3[4],
//     //                                     i);
//     // }
//     // std::cout << sum << " "
//     //           << test_points[3].x + test_points[3].y + test_points[3].z <<
//     //           '\n';

//     // sum = 0;
//     // for (int i = 0; i < 5; i++) {
//     //     // sum +=
//     //     //     getLinearBasisFunction11({0.25, 0.5, 0.2}, nodes_0[0],
//     //     //     nodes_0[1],
//     //     //                              nodes_0[2], nodes_0[3],
//     nodes_0[4],
//     //     i);

//     //     sum +=
//     //         result[indexes_0[i]] *
//     //         getLinearBasisFunction11({0.25, 0.5, 0.2}, nodes_0[0],
//     //         nodes_0[1],
//     //                                  nodes_0[2], nodes_0[3], nodes_0[4],
//     i);
//     // }
//     // std::cout << sum << " " << 0.25 + 0.5 + 0.2 << '\n';

//     finite_elements.push_back(
//                         {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 28});
//                     finite_elements.push_back(
//                         {0, 1, 2, 12, 13, 20, 21, 22, 8, 9, 16, 17, 28});
//                     finite_elements.push_back(
//                         {0, 3, 5, 12, 14, 20, 23, 25, 8, 10, 16, 18, 28});
//                     finite_elements.push_back(
//                         {5, 6, 7, 14, 15, 25, 26, 27, 10, 11, 18, 19, 28});
//                     finite_elements.push_back(
//                         {2, 4, 7, 13, 15, 22, 24, 27, 9, 11, 17, 19, 28});
//                     finite_elements.push_back(
//                         {20, 21, 22, 23, 24, 25, 26, 27, 16, 17, 18, 19,
//                         28});
//                     nodes.push_back({0, 0, 0});
//                     nodes.push_back({0.5, 0, 0});
//                     nodes.push_back({1, 0, 0});
//                     nodes.push_back({0, 0.5, 0});
//                     nodes.push_back({1, 0.5, 0});
//                     nodes.push_back({0, 1, 0});
//                     nodes.push_back({0.5, 1, 0});
//                     nodes.push_back({1, 1, 0});
//                     nodes.push_back({0.25, 0.25, 0.25});
//                     nodes.push_back({0.75, 0.25, 0.25});
//                     nodes.push_back({0.25, 0.75, 0.25});
//                     nodes.push_back({0.75, 0.75, 0.25});
//                     nodes.push_back({0, 0, 0.5});
//                     nodes.push_back({1, 0, 0.5});
//                     nodes.push_back({0, 1, 0.5});
//                     nodes.push_back({1, 1, 0.5});
//                     nodes.push_back({0.25, 0.25, 0.75});
//                     nodes.push_back({0.75, 0.25, 0.75});
//                     nodes.push_back({0.25, 0.75, 0.75});
//                     nodes.push_back({0.75, 0.75, 0.75});
//                     nodes.push_back({0, 0, 1});
//                     nodes.push_back({0.5, 0, 1});
//                     nodes.push_back({1, 0, 1});
//                     nodes.push_back({0, 0.5, 1});
//                     nodes.push_back({1, 0.5, 1});
//                     nodes.push_back({0, 1, 1});
//                     nodes.push_back({0.5, 1, 1});
//                     nodes.push_back({1, 1, 1});
//                     nodes.push_back({0.5, 0.5, 0.5});
   // auto integrateZ = [&](double x, double y) {
    //     return integrator.integrate(
    //         [&](double z) {
    //             double physical_x = 0.0;
    //             double physical_y = 0.0;
    //             double physical_z = 0.0;

    //             MathUtils::calculatePhysicalCoordinates(
    //                 {x, y, z}, node_0, node_1, node_2, node_3, vertex,
    //                 physical_x, physical_y, physical_z);

    //             return MathUtils::calculateF(
    //                        {physical_x, physical_y, physical_z}) *
    //                    MathUtils::getBasisFunction(
    //                        nodes, {x, y, z}, node_indexes,
    //                        basis_functions_type,
    //                        basis_functions_element_type,
    //                        index_of_basis_function, false);
    //         },
    //         0.0, 1.0);
    // };

    // auto integrateY = [&](double x) {
    //     return integrator.integrate([&](double y) { return integrateZ(x, y);
    //     },
    //                                 0.0, 1.0);
    // };

    // double result =
    //     integrator.integrate([&](double x) { return integrateY(x); },
    //     0.0, 1.0);


    // for (int i = 0; i < number_of_vertices_of_pyramid; i++) {
    //     for (int j = 0; j < number_of_vertices_of_pyramid; j++) {
    //         right_part_local[i] +=
    //             Integrator::integrateForMassMatrix(
    //                 nodes, finite_elements,
    //                 std::make_tuple(index_of_finite_element, i, j),
    //                 basis_functions_type, basis_functions_element_type) *
    //             MathUtils::calculateF(nodes[node_indexes.at(j)]);
    //     }
    // }

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

// Тестовая функция: x^4 + y^4 + z^4 (на эталонной пирамиде z = theta)
double testFunction(double xi, double eta, double theta) {
    return 1.0;
    // return std::pow(xi, 1) + std::pow(eta, 1) + std::pow(theta, 1);
}

double integrate(int n) {
    double dx = 1.0 / n;  // Шаг по x
    double dy = 1.0 / n;  // Шаг по y
    double dz;            // Шаг по z будет зависеть от x и y
    double result = 0.0;  // Итоговый интеграл

    for (int i = 0; i < n; ++i) {
        double x = i * dx + dx / 2.0;  // Центр прямоугольника по x
        for (int j = 0; j < n; ++j) {
            double y = j * dy + dy / 2.0;  // Центр прямоугольника по y
            double z_bound =
                1.0 -
                2.0 * std::max(std::abs(x - 0.5),
                               std::abs(y - 0.5));  // Верхняя граница для z
            dz = z_bound / n;                       // Делим на n для z

            for (int k = 0; k < n; ++k) {
                double z = k * dz + dz / 2.0;  // Центр прямоугольника по z
                result += testFunction(x, y, z) * dx * dy * dz;
            }
        }
    }

    return result;
}

int main() {
    // Интеграл тестовой функции
    double integral = integrate(100);
    std::cout << "Integral of x^4 + y^4 + z^4: " << integral << std::endl;
    return 0;
}
