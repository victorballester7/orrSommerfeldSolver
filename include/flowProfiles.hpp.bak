#ifndef FLOW_PROFILES_HPP
#define FLOW_PROFILES_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Eigenvalues>
#include <stdexcept>
#include <string>
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
  double a, b;
  double r;
  bool useRationalMapping;

  void readFromFile(const std::string &filename, std::vector<double> &xdata,
                    uint colX, std::vector<double> &data, uint colY,
                    uint numSkipHeaderLines) const;

  double interpolate(double z, const std::vector<double> &xdata,
                     const std::vector<double> &ydata) const;
  double interpolateSpline(double z, const std::vector<double> &xdata,
                           const std::vector<double> &ydata,
                           const std::vector<double> &y2data) const;
  std::vector<double>
  buildNaturalSplineSecondDerivatives(const std::vector<double> &xdata,
                                      const std::vector<double> &ydata) const;

  std::vector<double> diffData(const std::vector<double> &xdata,
                               const std::vector<double> &ydata) const;

public:
  std::vector<double> x_data;
  std::vector<double> u_data;
  std::vector<double> du_data;
  std::vector<double> d2u_data;
  std::vector<double> u_spline2;
  std::vector<double> d2u_spline2;

  CustomU(std::string filename, uint colX, uint colY, uint numSkipHeaderLines,
          bool useRationalMapping_ = false)
      : useRationalMapping(useRationalMapping_) {
    readFromFile(filename, x_data, colX, u_data, colY, numSkipHeaderLines);
    if (x_data.size() < 3 || u_data.size() < 3) {
      throw std::invalid_argument(
          "CustomU profile must contain at least 3 points.");
    }
    if (x_data.size() != u_data.size()) {
      throw std::invalid_argument(
          "CustomU profile x and U columns have inconsistent lengths.");
    }
    if (!std::is_sorted(x_data.begin(), x_data.end())) {
      throw std::invalid_argument(
          "CustomU x_data must be sorted in ascending order.");
    }
    if (useRationalMapping) {
      // only use the portion of the data near the wall and forget out the far
      // field
      const double maxY = x_data.front() + 25.0;
      auto upper = std::upper_bound(x_data.begin(), x_data.end(), maxY);
      const size_t keep =
          static_cast<size_t>(std::distance(x_data.begin(), upper));
      if (keep >= 3 && keep < x_data.size()) {
        x_data.resize(keep);
        u_data.resize(keep);
      }
    }
    a = x_data[0];
    b = x_data.back();
    if (b <= a) {
      throw std::invalid_argument(
          "CustomU profile requires strictly increasing x range.");
    }
    const double span = b - a;
    if (useRationalMapping) {
      double R = -1.0;
      // Place the clustering point near the edge of the shear layer.
      const double Uinf = u_data.back();
      const double target = 0.995 * Uinf;
      uint idx = static_cast<uint>(u_data.size() - 1);
      for (uint i = 0; i < static_cast<uint>(u_data.size()); ++i) {
        if (u_data[i] >= target) {
          idx = i;
          break;
        }
      }
      R = x_data[idx] - a;
      const double minR = 0.02 * span;
      const double maxR = 0.45 * span;
      R = std::clamp(R, minR, maxR);
      r = a + R;
    }

    du_data = diffData(x_data, u_data);
    d2u_data = diffData(x_data, du_data);
    refreshInterpolationData();
  }

  double mapToStandardRegion(double y) const override {
    if (useRationalMapping) {
      const double B = b - a;
      const double R = r - a;
      const double t = y - a;
      return B * (t - R) / (R * B + t * (B - 2.0 * R));
    } else {
      return 2.0 * (y - a) / (b - a) - 1.0;
    }
  }

  double mapToPhysicalRegion(double xi) const override {
    if (useRationalMapping) {
      const double B = b - a;
      const double R = r - a;
      return a + R * B * (1.0 + xi) / (B - xi * (B - 2.0 * R));
    } else {
      return (b + a) / 2.0 + (b - a) / 2.0 * xi;
    }
  }

  std::vector<double>
  getJacobian(std::vector<double> gaussPoints) const override {
    if (useRationalMapping) {
      const double B = b - a;
      const double R = r - a;
      const double denomConst = 2.0 * B * R * (B - R);
      std::vector<double> jacobian(gaussPoints.size(), 0.0);

      for (uint i = 0; i < gaussPoints.size(); i++) {
        const double den = B - gaussPoints[i] * (B - 2.0 * R);
        jacobian[i] = den * den / denomConst;
      }
      return jacobian;
    } else {
      // Compute the Jacobian of the mapping [a, b] -> [-1, 1]
      double jac = 2.0 / (b - a);
      std::vector<double> jacobian(gaussPoints.size(), jac);
      return jacobian;
    }
  }

  std::vector<double>
  getDiffJacobian(std::vector<double> gaussPoints) const override {
    // Return d/dy(dxi/dy), needed by chain rule in d2N.
    std::vector<double> djacobian(gaussPoints.size(), 0.0);
    if (useRationalMapping) {
      const double B = b - a;
      const double R = r - a;
      const double denomConst = 2.0 * B * B * R * R * (B - R) * (B - R);
      for (uint i = 0; i < gaussPoints.size(); i++) {
        const double den = B - gaussPoints[i] * (B - 2.0 * R);
        djacobian[i] = -std::pow(den, 3) * (B - 2.0 * R) / denomConst;
      }
      return djacobian;
    } else {
      // The Jacobian is constant, so its derivative is zero.
      return djacobian;
    }
  }

  double getU(double z) const override {
    return interpolateSpline(z, x_data, u_data, u_spline2);
  }

  double getd2U(double z) const override {
    return interpolateSpline(z, x_data, d2u_data, d2u_spline2);
  }

  void refreshInterpolationData() {
    u_spline2 = buildNaturalSplineSecondDerivatives(x_data, u_data);
    d2u_spline2 = buildNaturalSplineSecondDerivatives(x_data, d2u_data);
  }
};

#endif // FLOW_PROFILES_HPP
