#include <autodiff/forward/dual.hpp>

using namespace autodiff;

dual getDualLinearBasisFunctions(int basis_function_index, dual xi, dual eta,
                                 dual theta) {
    dual xi_0 = 1 - xi;
    dual xi_1 = xi;

    dual eta_0 = 1 - eta;
    dual eta_1 = eta;

    dual theta_0 = 1 - theta;
    dual theta_1 = theta;

    switch (basis_function_index) {
        case 0:
            return xi_0 * eta_0 * theta_0;
            break;
        case 1:
            return xi_1 * eta_0 * theta_0;
            break;
        case 2:
            return xi_0 * eta_1 * theta_0;
            break;
        case 3:
            return xi_1 * eta_1 * theta_0;
            break;
        case 4:
            return theta_1;
            break;
        default:
            return 0.0;
            break;
    }

    return 0.0;
}

dual getDualQuadraticBasisFunctions(int basis_function_index, dual xi, dual eta,
                                    dual theta) {
    dual xi_0 = (1 - 2 * xi) * (1 - xi);
    dual xi_1 = 4 * xi * (1 - xi);
    dual xi_2 = -xi * (1 - 2 * xi);

    dual eta_0 = (1 - 2 * eta) * (1 - eta);
    dual eta_1 = 4 * eta * (1 - eta);
    dual eta_2 = -eta * (1 - 2 * eta);

    dual theta_0 = (1 - 2 * theta) * (1 - theta);
    dual theta_1 = 4 * theta * (1 - theta);
    dual theta_2 = -theta * (1 - 2 * theta);

    switch (basis_function_index) {
        case 0:
            return xi_0 * eta_0 * theta_0;
            break;
        case 1:
            return xi_1 * eta_0 * theta_0;
            break;
        case 2:
            return xi_2 * eta_0 * theta_0;
            break;
        case 3:
            return xi_0 * eta_1 * theta_0;
            break;
        case 4:
            return xi_1 * eta_1 * theta_0;
            break;
        case 5:
            return xi_2 * eta_1 * theta_0;
            break;
        case 6:
            return xi_0 * eta_2 * theta_0;
            break;
        case 7:
            return xi_1 * eta_2 * theta_0;
            break;
        case 8:
            return xi_2 * eta_2 * theta_0;
            break;
        case 9:
            return 4 * (0.75 - xi) * (0.75 - eta) * theta_1;
            break;
        case 10:
            return -4 * (0.25 - xi) * (0.75 - eta) * theta_1;
            break;
        case 11:
            return -4 * (0.75 - xi) * (0.25 - eta) * theta_1;
            break;
        case 12:
            return 4 * (0.25 - xi) * (0.25 - eta) * theta_1;
            break;
        case 13:
            return theta_2;
            break;
        default:
            return 0.0;
            break;
    }

    return 0.0;
}

dual getDualCubicBasisFunctions(int basis_function_index, dual xi, dual eta,
                                dual theta) {
    dual xi_0 = 1.0 / 2.0 * (1 - 3 * xi) * (2 - 3 * xi) * (1 - xi);
    dual xi_1 = 9.0 / 2.0 * xi * (2 - 3 * xi) * (1 - xi);
    dual xi_2 = -9.0 / 2.0 * xi * (1 - 3 * xi) * (1 - xi);
    dual xi_3 = 1.0 / 2.0 * xi * (1 - 3 * xi) * (2 - 3 * xi);

    dual eta_0 = 1.0 / 2.0 * (1 - 3 * eta) * (2 - 3 * eta) * (1 - eta);
    dual eta_1 = 9.0 / 2.0 * eta * (2 - 3 * eta) * (1 - eta);
    dual eta_2 = -9.0 / 2.0 * eta * (1 - 3 * eta) * (1 - eta);
    dual eta_3 = 1.0 / 2.0 * eta * (1 - 3 * eta) * (2 - 3 * eta);

    dual theta_0 = 1.0 / 2.0 * (1 - 3 * theta) * (2 - 3 * theta) * (1 - theta);
    dual theta_1 = 9.0 / 2.0 * theta * (2 - 3 * theta) * (1 - theta);
    dual theta_2 = -9.0 / 2.0 * theta * (1 - 3 * theta) * (1 - theta);
    dual theta_3 = 1.0 / 2.0 * theta * (1 - 3 * theta) * (2 - 3 * theta);

    switch (basis_function_index) {
        case 0:
            return xi_0 * eta_0 * theta_0;
            break;
        case 1:
            return xi_1 * eta_0 * theta_0;
            break;
        case 2:
            return xi_2 * eta_0 * theta_0;
            break;
        case 3:
            return xi_3 * eta_0 * theta_0;
            break;
        case 4:
            return xi_0 * eta_1 * theta_0;
            break;
        case 5:
            return xi_1 * eta_1 * theta_0;
            break;
        case 6:
            return xi_2 * eta_1 * theta_0;
            break;
        case 7:
            return xi_3 * eta_1 * theta_0;
            break;
        case 8:
            return xi_0 * eta_2 * theta_0;
            break;
        case 9:
            return xi_1 * eta_2 * theta_0;
            break;
        case 10:
            return xi_2 * eta_2 * theta_0;
            break;
        case 11:
            return xi_3 * eta_2 * theta_0;
            break;
        case 12:
            return xi_0 * eta_3 * theta_0;
            break;
        case 13:
            return xi_1 * eta_3 * theta_0;
            break;
        case 14:
            return xi_2 * eta_3 * theta_0;
            break;
        case 15:
            return xi_3 * eta_3 * theta_0;
            break;
        case 16:
            return 1.0 / 16.0 * (5 - 6 * xi) * (5 - 6 * eta) * theta_1;
            break;
        case 17:
            return -1.0 / 16.0 * (1 - 6 * xi) * (5 - 6 * eta) * theta_1;
            break;
        case 18:
            return -1.0 / 16.0 * (5 - 6 * xi) * (1 - 6 * eta) * theta_1;
            break;
        case 19:
            return 1.0 / 16.0 * (1 - 6 * xi) * (1 - 6 * eta) * theta_1;
            break;
        case 20:
            return (2 - 3 * xi) * (2 - 3 * eta) * theta_2;
            break;
        case 21:
            return -(1 - 3 * xi) * (2 - 3 * eta) * theta_2;
            break;
        case 22:
            return -(2 - 3 * xi) * (1 - 3 * eta) * theta_2;
            break;
        case 23:
            return (1 - 3 * xi) * (1 - 3 * eta) * theta_2;
            break;
        case 24:
            return theta_3;
            break;
        default:
            return 0.0;
            break;
    }

    return 0.0;
}

dual getDualHierarhicalQuadraticBasisFunctions(int basis_function_index,
                                               dual xi, dual eta, dual theta) {
    dual xi_0 = 1 - xi;
    dual xi_1 = xi;
    dual xi_2 = 1.0 / 8.0 * ((2 * xi - 1) * (2 * xi - 1) - 1);

    dual eta_0 = 1 - eta;
    dual eta_1 = eta;
    dual eta_2 = 1.0 / 8.0 * ((2 * eta - 1) * (2 * eta - 1) - 1);

    dual theta_0 = 1 - theta;
    dual theta_1 = theta;
    dual theta_2 = 1.0 / 8.0 * ((2 * theta - 1) * (2 * theta - 1) - 1);

    switch (basis_function_index) {
        case 0:
            return xi_0 * eta_0 * theta_0;
            break;
        case 1:
            return xi_2 * eta_0 * theta_0;
            break;
        case 2:
            return xi_1 * eta_0 * theta_0;
            break;
        case 3:
            return xi_0 * eta_2 * theta_0;
            break;
        case 4:
            return xi_2 * eta_2 * theta_0;
            break;
        case 5:
            return xi_1 * eta_2 * theta_0;
            break;
        case 6:
            return xi_0 * eta_1 * theta_0;
            break;
        case 7:
            return xi_2 * eta_1 * theta_0;
            break;
        case 8:
            return xi_1 * eta_1 * theta_0;
            break;
        case 9:
            return xi_0 * eta_0 * theta_2;
            break;
        case 10:
            return xi_1 * eta_0 * theta_2;
            break;
        case 11:
            return xi_0 * eta_1 * theta_2;
            break;
        case 12:
            return xi_1 * eta_1 * theta_2;
            break;
        case 13:
            return theta_1;
            break;
        default:
            return 0.0;
            break;
    }

    return 0.0;
}

dual getDualHierarhicalCubicBasisFunctions(int basis_function_index, dual xi,
                                           dual eta, dual theta) {
    return 0.0;
}

dual getDualBasisFunction(BASIS_TYPE basis_functions_type,
                          BASIS_ELEMENT_TYPE basis_functions_elements_type,
                          int basis_function_index, dual xi, dual eta,
                          dual theta) {
    switch (basis_functions_type) {
        case BASIS_TYPE::Lagrange:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Linear:
                    return getDualLinearBasisFunctions(basis_function_index, xi,
                                                       eta, theta);
                    break;
                case BASIS_ELEMENT_TYPE::Quadratic:
                    return getDualQuadraticBasisFunctions(basis_function_index,
                                                          xi, eta, theta);
                    break;
                case BASIS_ELEMENT_TYPE::Cubic:
                    return getDualCubicBasisFunctions(basis_function_index, xi,
                                                      eta, theta);
                    break;
                default:
                    return 0.0;
                    break;
            }
            break;
        case BASIS_TYPE::Hierarhical:
            switch (basis_functions_elements_type) {
                case BASIS_ELEMENT_TYPE::Quadratic:
                    return getDualHierarhicalQuadraticBasisFunctions(
                        basis_function_index, xi, eta, theta);
                    break;
                case BASIS_ELEMENT_TYPE::Cubic:
                    return getDualHierarhicalCubicBasisFunctions(
                        basis_function_index, xi, eta, theta);
                    break;
                default:
                    return 0.0;
                    break;
            }
            break;
        default:
            return 0.0;
            break;
    }

    return 0.0;
}

void computeDerivative(BASIS_TYPE basis_functions_type,
                       BASIS_ELEMENT_TYPE basis_functions_elements_type,
                       int basis_function_index, dual ksi, dual nu, dual tetta,
                       double& d_psi_d_xi, double& d_psi_d_eta,
                       double& d_psi_d_theta) {
    auto result = derivatives(
        [&](dual xi, dual eta, dual theta) {
            return getDualBasisFunction(basis_functions_type,
                                        basis_functions_elements_type,
                                        basis_function_index, xi, eta, theta);
        },
        wrt(ksi), at(ksi, nu, tetta));
    d_psi_d_xi = std::get<1>(result);

    result = derivatives(
        [&](dual xi, dual eta, dual theta) {
            return getDualBasisFunction(basis_functions_type,
                                        basis_functions_elements_type,
                                        basis_function_index, xi, eta, theta);
        },
        wrt(nu), at(ksi, nu, tetta));
    d_psi_d_eta = std::get<1>(result);

    result = derivatives(
        [&](dual xi, dual eta, dual theta) {
            return getDualBasisFunction(basis_functions_type,
                                        basis_functions_elements_type,
                                        basis_function_index, xi, eta, theta);
        },
        wrt(tetta), at(ksi, nu, tetta));
    d_psi_d_theta = std::get<1>(result);
}
