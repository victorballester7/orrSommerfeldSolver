#ifndef FLOW_PROFILES_HPP
#define FLOW_PROFILES_HPP

#include <cstddef>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <vector>

class Uprofile {
public:
  // Flow profile: U(z) for plane Poiseuille flow
  virtual double getU(double z) const = 0;

  // Flow profile second derivative: U''(z)
  virtual double getd2U(double z) const = 0;

  virtual double mapToStandardRegion(double y) const {
    // Map the flow profile from the physical region [a, b] to the standard
    // region [-1, 1]
    return y;
  }

  // Map the flow profile from the standard region [-1, 1] to the physical
  // region [a, b]
  virtual double mapToPhysicalRegion(double xi) const {
    // Map the flow profile from the standard region [-1, 1] to the physical
    // region [a, b]
    return xi;
  }

  // Jacobian of the mapping [a, b] -> [-1, 1]
  virtual std::vector<double>
  getJacobian(std::vector<double> gaussPoints) const {
    // Jacobian of the mapping [-1, 1] -> [-1, 1]
    std::vector<double> jacobian(gaussPoints.size(), 1.0);
    return jacobian;
  }

  // Derivative of the Jacobian of the mapping [a, b] -> [-1, 1]
  virtual std::vector<double>
  getDiffJacobian(std::vector<double> gaussPoints) const {
    // Derivative of the Jacobian of the mapping [-1, 1] -> [-1, 1]
    std::vector<double> djacobian(gaussPoints.size(), 0.0);
    return djacobian;
  }
};

// Plane Poiseuille flow profile
class Poiseuille : public Uprofile {
public:
  Poiseuille() = default;

  double getU(double z) const override {
    return 1.0 - z * z; // Plane Poiseuille flow
  }

  double getd2U(__attribute__((unused)) double z) const override {
    return -2.0;
  }

  double mapToStandardRegion(double y) const override {
    return Uprofile::mapToStandardRegion(y);
  }

  double mapToPhysicalRegion(double xi) const override {
    return Uprofile::mapToPhysicalRegion(xi);
  }

  std::vector<double>
  getJacobian(std::vector<double> gaussPoints) const override {
    return Uprofile::getJacobian(gaussPoints);
  }

  std::vector<double>
  getDiffJacobian(std::vector<double> gaussPoints) const override {
    return Uprofile::getDiffJacobian(gaussPoints);
  }
};

// Plane Couette flow profile
class Couette : public Uprofile {
public:
  Couette() = default;

  double getU(double z) const override {
    return z; // Plane Couette flow
  }

  double getd2U(__attribute__((unused)) double z) const override { return 0.0; }

  double mapToStandardRegion(double y) const override {
    return Uprofile::mapToStandardRegion(y);
  }

  double mapToPhysicalRegion(double xi) const override {
    return Uprofile::mapToPhysicalRegion(xi);
  }

  std::vector<double>
  getJacobian(std::vector<double> gaussPoints) const override {
    return Uprofile::getJacobian(gaussPoints);
  }

  std::vector<double>
  getDiffJacobian(std::vector<double> gaussPoints) const override {
    return Uprofile::getDiffJacobian(gaussPoints);
  }
};

// Blasius flow profile
class CustomU : public Uprofile {
private:
  size_t len_data;
  double a, b;
  double r;

  void readFromFile(const std::string &filename, std::vector<double> &xdata,
                    uint colX, std::vector<double> &data, uint colY,
                    uint numSkipHeaderLines) const;

  double interpolate(double z, const std::vector<double> &xdata,
                     const std::vector<double> &ydata) const;

  std::vector<double> diffData(const std::vector<double> &xdata,
                               const std::vector<double> &ydata) const;

public:
  std::vector<double> x_data;
  std::vector<double> u_data;
  std::vector<double> du_data;
  std::vector<double> d2u_data;

  CustomU(std::string filename, uint colX, uint colY, uint numSkipHeaderLines) {
    readFromFile(filename, x_data, colX, u_data, colY, numSkipHeaderLines);
    a = x_data[0];
    b = x_data.back();
    r = a + 60; // half of the points are clustered between a and r
    len_data = u_data.size();
    du_data = diffData(x_data, u_data);
    d2u_data = diffData(x_data, du_data);
  }

  double mapToStandardRegion(double y) const override {
    return 2.0 * (y - a) / (b - a) - 1.0;
    // return b * (y - r) / (r * b + y * (b - 2 * r));
  }

  double mapToPhysicalRegion(double xi) const override {
    return (b + a) / 2.0 + (b - a) / 2.0 * xi;
    // return r * b * (1 + xi) / (b - xi * (b - 2 * r));
  }

  std::vector<double>
  getJacobian(std::vector<double> gaussPoints) const override {
    // Compute the Jacobian of the mapping [a, b] -> [-1, 1]
    double jac = 2.0 / (b - a);
    std::vector<double> jacobian(gaussPoints.size(), jac);
    

    // double numerator = 2. * b * r * (b - r);
    // std::vector<double> jacobian(gaussPoints.size(), 0.0);

    // for (uint i = 0; i < gaussPoints.size(); i++) {
    //   // jacobian[i] = 2.0 / (b - a);
    //   jacobian[i] =
    //   numerator / std::pow(r * b + gaussPoints[i] * (b - 2 * r), 2);
    // }
    return jacobian;
  }

  std::vector<double>
  getDiffJacobian(std::vector<double> gaussPoints) const override {
    // Compute the derivative of the Jacobian of the mapping [a, b] -> [-1, 1]
    std::vector<double> djacobian(gaussPoints.size(), 0.0);

    // double numerator = -4. * b * r * (b - r) * (b - 2 * r);
    // std::vector<double> jacobian(gaussPoints.size(), 0.0);

    // for (uint i = 0; i < gaussPoints.size(); i++) {
    //   // jacobian[i] = 0.0;
    //   jacobian[i] =
    //       numerator / std::pow(r * b + gaussPoints[i] * (b - 2 * r), 3);
    // }
    return djacobian;
  }

  double getU(double z) const override {
    return interpolate(z, x_data, u_data);
  }

  double getd2U(double z) const override {
    return interpolate(z, x_data, d2u_data);
  }
};

#endif // FLOW_PROFILES_HPP
