#ifndef FAST_CIRCLE_FIT_H
#define FAST_CIRCLE_FIT_H

#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <cmath>
#include <iostream>

// Simple point structure (x,y)
struct Point2D_xy {
    double x;
    double y;
};

// Simple error structure: angular error sigma_phi
struct HitError {
    double sigma_phi; 
};

/*
 * Generic version of the CMS FastCircleFit.
 * Fits a circle to a set of weighted 2D hits using the method
 * described by Strandlie & Frühwirth (NIM A488 (2002)).
*/


class FastCircleFit {

public:

    FastCircleFit(const std::vector<Point2D_xy>& points)// , const std::vector<HitError>& errors = {} )
    {
        size_t N = points.size();
        std::vector<double> x(N), y(N), z(N), w(N);
        // calculate(points, errors, x, y, z, w);
        calculate(points, x, y, z, w);
    }


    Point2D_xy closestPointTo(const Point2D_xy& hit) const
    {
        // vector from circle center to hit
        double ux = hit.x - x0_;
        double uy = hit.y - y0_;

        double r = std::sqrt(ux*ux + uy*uy);

        // if the hit is exactly at the center (degenerate)
        if (r == 0) {
            // return arbitrary point on the circle
            return Point2D_xy{x0_ + rho_, y0_};
        }

        double scale = rho_ / r;

        // closest point on circle
        return Point2D_xy{
            x0_ + scale * ux,
            y0_ + scale * uy
        };
    }

    Point2D_xy tangentAtPCA(const Point2D_xy& pca, const Point2D_xy& nextHit) const
    {

        Point2D_xy tangent;

        double dx = pca.x - x0_;
        double dy = pca.y - y0_;
        double dx_next = nextHit.x - x0_;
        double dy_next = nextHit.y - y0_;

        // cross product Z component
        double crossZ = dx*dy_next - dy*dx_next;
        if (crossZ < 0) {

            tangent.x =  dy;
            tangent.y = -dx;
        } else {
          
            tangent.x = -dy;
            tangent.y =  dx;
        }

        double norm = std::sqrt(tangent.x * tangent.x + tangent.y * tangent.y);
        if (norm > 0) {
            tangent.x /= norm;
            tangent.y /= norm;
        }

        return tangent;
    }


    double x0() const { return x0_; }
    double y0() const { return y0_; }
    double rho() const { return rho_; }
    double chi2() const { return chi2_; }

private:

    template <typename C>
    // void calculate(const std::vector<Point2D_xy>& pts,
    //                const std::vector<HitError>& errs,
    //                C& x, C& y, C& z, C& weight)
    void calculate(const std::vector<Point2D_xy>& pts, C& x, C& y, C& z, C& weight)
    {
        const size_t N = pts.size();

        // Transform input points
        for (size_t i = 0; i < N; ++i) {
            x[i] = pts[i].x;
            y[i] = pts[i].y;
            z[i] = pts[i].x * pts[i].x + pts[i].y * pts[i].y;  // artificial third coordinate

            // Weight = 1 / (r^2 * sigma_phi^2)
            double r2 = x[i] * x[i] + y[i] * y[i];
            // float sigma_phi2 = errs[i].sigma_phi * errs[i].sigma_phi;
            // weight[i] = 1.0f / (r2 * sigma_phi2);
            weight[i] = 1.0f / (r2);
        }


        double totalWeight = std::accumulate(weight.begin(), weight.end(), 0.0f);
        double invTotWeight = 1.0f / totalWeight;

        // Weighted mean in (x, y, z)
        Eigen::Vector3d mean = Eigen::Vector3d::Zero();
        for (size_t i = 0; i < N; ++i) {
            mean[0] += weight[i] * x[i];
            mean[1] += weight[i] * y[i];
            mean[2] += weight[i] * z[i];
        }
        mean *= invTotWeight;

        // Build covariance matrix
        Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
        for (size_t i = 0; i < N; ++i) {
            double dx = x[i] - mean[0];
            double dy = y[i] - mean[1];
            double dz = z[i] - mean[2];
            double w  = weight[i];

            A(0,0) += w * dx * dx;
            A(1,0) += w * dy * dx;
            A(1,1) += w * dy * dy;
            A(2,0) += w * dz * dx;
            A(2,1) += w * dz * dy;
            A(2,2) += w * dz * dz;
        }
        A *= 1.0f / double(N);

        // Eigenvalue decomposition
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(A);
        Eigen::Vector3d n = eigen.eigenvectors().col(0);  // smallest eigenvalue

        // Plane: n•(x,y,z) + c = 0
        double c = - (n[0]*mean[0] + n[1]*mean[1] + n[2]*mean[2]);

        // Extract circle parameters
        double tmp = 0.5 / n[2];
        x0_ = -n[0] * tmp;
        y0_ = -n[1] * tmp;

        double rho2 = (1 - n[2]*n[2] - 4*c*n[2]) * tmp * tmp;
        rho_ = std::sqrt(rho2);

        // Compute chi2-like sum
        double S = 0;
        for (size_t i = 0; i < N; ++i) {
            double res = c + n[0]*x[i] + n[1]*y[i] + n[2]*z[i];
            S += res*res * weight[i];
        }
        chi2_ = S;
    }

    double x0_{0};
    double y0_{0};
    double rho_{0};
    double chi2_{0};
};

#endif // FAST_CIRCLE_FIT_H
