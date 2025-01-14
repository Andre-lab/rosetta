// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       core/scoring/shape/zernike/Zernike2Ddescriptor.cc
/// @brief      2D Zernike transform class
/// @details    ... TO FILL
/// @author     Ingemar Andre

#include <core/scoring/shape/zernike/Zernike2Ddescriptor.hh>
#include <basic/Tracer.hh>
#include <iostream>
#include <math.h>
#include <chrono>

namespace core {
namespace scoring {
namespace shape {
namespace zernike {

static basic::Tracer TR( "core.scoring.shape.zernike" );

//using namespace std::complex_literals;

double Zernike2Ddescriptor::get_r(int x,int y, int N){
    int center=N/2;
    double r2 = pow(x - center,2) + pow(y - center,2);
    double r = sqrt(r2);
    return r;
}

unsigned long int
Zernike2Ddescriptor::factorial(int n) {

    unsigned long long factorial = 1;
    for(int i = 1; i <=n; ++i) {
        factorial *= i;
    }
    return factorial;
}

int Zernike2Ddescriptor::count_nonzero(std::vector<std::vector<double>> & rho) {
    int num(0);
    for (int xi=0; xi< int(rho.size()); xi++) {
        for (int yi=0; yi < int(rho[xi].size()); yi++) {
            if (rho[xi][yi] > 0) ++num;
        }
    }
    return num;
}

std::vector< std::vector<int> >
Zernike2Ddescriptor::get_nm(int n) {
    std::vector<std::vector<int> > nm;
    for (int ni = 0; ni <=n; ni++) {
        std::vector<int> mv;
//        for (int mi = -ni; mi <= ni; mi++) {
        for (int mi = 0; mi <= ni; mi++) {
            int absm = abs(mi);
            if ((ni - absm) % 2 == 0) {
                mv.push_back(mi);
            }
        }
        nm.push_back(mv);
    }
    return nm;
}

std::vector<std::vector<double>>
Zernike2Ddescriptor::rad_poly(std::vector<std::vector<double>> &rho,int n,int m) {
    int N = rho.size();
    int absm = abs(m);
    std::vector<std::vector<double>> sum_r;
    for (int i=0; i<N; i++) sum_r.push_back(std::vector<double>(N));

    for (int i=0; i<N; i++) {
        for (int j = 0; j < N; j++) {
            for (int s = 0; s <= int((n - absm) / 2); s++) {
                double val = pow(-1, s) * factorial(n - s) /
                             (factorial(s) * factorial((n + absm) / 2 - s) * factorial((n - absm) / 2 - s));
                sum_r[i][j] = sum_r[i][j] + val * pow(rho[i][j], n - 2 * s);
            }
        }
    }
    return sum_r;
}

Zernike2Ddescriptor::Zernike2Ddescriptor(){
	initialized_ = false;
}


Zernike2Ddescriptor::Zernike2Ddescriptor(std::vector<std::vector<int>> grid, int grid_size, int order) {
    grid_ = grid;
    dim_ = grid_size;
    order_ = order;
    std::vector<std::vector<double>> rho;
    std::vector<std::vector<double>> rad;
    std::vector<std::vector<double>> theta;
//    std::vector<std::vector<int>> grid2;
    std::vector<std::vector< std::complex<double> >> Zprod;
    for (int i=0; i<dim_; i++) {
//        grid2.push_back(std::vector<int>(dim_));
        rho.push_back(std::vector<double>(dim_));
        rad.push_back(std::vector<double>(dim_));
        theta.push_back(std::vector<double>(dim_));
        Zprod.push_back(std::vector< std::complex<double> >(dim_));
    }
/*
      double radius = 40;
      double a=20;
      double b=10;
      for (int i=0; i<dim_; i++) {
        for (int j=0; j<dim_; j++) {
            if ( pow(i-dim_/2,2)/pow(a,2) + pow(j-dim_/2,2)/pow(b,2) < 1 ) {
//            if (get_r(i,j,dim_) < radius) {
                grid2[i][j] = 1;
            }
        }
    }
*/
    std::vector<int> X,Y;
    for (int i=0; i<dim_; i++) {
        X.push_back(i);
        Y.push_back(i);
    }

    for (int xi=0; xi<dim_; xi++) {
        for (int yi = 0; yi < dim_; yi++) {
            rho[xi][yi] = sqrt(pow(2*xi - dim_ + 1,2) + pow(2*yi - dim_ + 1,2))/dim_;
            if (rho[xi][yi] > 1){
                rho[xi][yi] = 0;
            }
        }
    }

    for (int xi=0; xi<dim_; xi++) {
        for (int yi = 0; yi < dim_; yi++) {
            theta[xi][yi] = atan2( (dim_-1-2*yi),(2*xi-dim_+1));
        }
    }

    auto start = std::chrono::high_resolution_clock::now();
    invariants_.clear();

    std::vector<std::vector<int>> nm = get_nm(order);
		if ( rad_poly_store_.size() == 0 ) {
//			std::cout << "Calling rad update..." << "order..." << order << std::endl;
			for (int ni=0; ni< int(nm.size());ni++ ) {
        	for (int j=0; j< int(nm[ni].size());j++ ) {
            int n = ni;
            int m = nm[ni][j];
						rad = rad_poly(rho,n,m);
						rad_poly_store_.push_back(rad);
				}
		}
	}
    int counter(0);
    for (int ni=0; ni< int(nm.size());ni++ ) {
        for (int j=0; j< int(nm[ni].size());j++ ) {
            int n = ni;
            int m = nm[ni][j];
						rad = rad_poly_store_[counter];
            counter++;
//            rad = rad_poly(rho,n,m);

            std::complex<double> Zsum(0);
            for (int xi=0; xi<dim_; xi++) {
                for (int yi = 0; yi < dim_; yi++) {
		                std::complex<double> imval(0, -theta[xi][yi] * double(m) );
//                    Zprod[xi][yi] = grid2[xi][yi] * rad[xi][yi] * exp(imval);
                    Zprod[xi][yi] = grid[xi][yi] * rad[xi][yi] * exp(imval);
                    Zsum += Zprod[xi][yi];
//                       std::cout << xi << " " << yi << " " << grid[xi][yi] << " " << rad[xi][yi] << " " << theta[xi][yi] << " " << Zprod[xi][yi] <<  std::endl;

                }
            }
            double norm_factor = double(n+1)/(count_nonzero(rho) + 1 );
            Zsum=Zsum*norm_factor;
            double Znorm = std::norm(Zsum);
            invariants_.push_back(std::sqrt(Znorm));
//            std::cout << Zsum << std::endl;
        }
    }
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
//    std::cout << "Time per transform: " << elapsed.count()  << std::endl;

  }

Zernike2Ddescriptor::Zernike2Ddescriptor(int grid_size, int order) {
//		std::cout << "Initializing Zernike2Ddescriptor..." << std::endl;
    dim_ = grid_size;
    order_ = order;
		setup();
    initialized_ = true;
}

bool Zernike2Ddescriptor::initialized() {

	return initialized_;

}

void
Zernike2Ddescriptor::setup() {

    nm_ = get_nm(order_);
//    std::vector<std::vector<int>> grid2;
    std::vector<std::vector< std::complex<double> >> Zprod;
    for (int i=0; i<dim_; i++) {
//        grid2.push_back(std::vector<int>(dim_));
        rho_.push_back(std::vector<double>(dim_));
//        rad_.push_back(std::vector<double>(dim_));
        theta_.push_back(std::vector<double>(dim_));
        Zprod.push_back(std::vector< std::complex<double> >(dim_));
    }
/*
      double radius = 40;
      double a=20;
      double b=10;
      for (int i=0; i<dim_; i++) {
        for (int j=0; j<dim_; j++) {
            if ( pow(i-dim_/2,2)/pow(a,2) + pow(j-dim_/2,2)/pow(b,2) < 1 ) {
//            if (get_r(i,j,dim_) < radius) {
                grid2[i][j] = 1;
            }
        }
    }
*/
    for (int xi=0; xi<dim_; xi++) {
        for (int yi = 0; yi < dim_; yi++) {
            rho_[xi][yi] = sqrt(pow(2*xi - dim_ + 1,2) + pow(2*yi - dim_ + 1,2))/dim_;
            if (rho_[xi][yi] > 1){
                rho_[xi][yi] = 0;
            }
        }
    }

    for (int xi=0; xi<dim_; xi++) {
        for (int yi = 0; yi < dim_; yi++) {
            theta_[xi][yi] = atan2( (dim_-1-2*yi),(2*xi-dim_+1));
        }
    }

//    auto start = std::chrono::high_resolution_clock::now();

		if ( rad_poly_store_.size() == 0 ) {
//			std::cout << "Calling rad update..." << "...with order..." << int(nm_.size()) << std::endl;
			for (int ni=0; ni< int(nm_.size());ni++ ) {
        	for (int j=0; j< int(nm_[ni].size());j++ ) {
          //  int n = ni;
            int m = nm_[ni][j];
						std::vector<std::vector<double>> rad = rad_poly(rho_,ni,m);
						rad_poly_store_.push_back(rad);
				}
		}
	}
}

void
Zernike2Ddescriptor::Transform( std::vector<std::vector<int>> grid ) {
//    auto start = std::chrono::high_resolution_clock::now();
    std::vector<std::vector< std::complex<double> >> Zprod;
    for (int i=0; i<dim_; i++) {
        Zprod.push_back(std::vector< std::complex<double> >(dim_));
    }

	std::vector<int> X,Y;
    for (int i=0; i<dim_; i++) {
        X.push_back(i);
        Y.push_back(i);
    }

    invariants_.clear();

	int counter(0);
    for (int ni=0; ni< int(nm_.size());ni++ ) {
        for (int j=0; j< int(nm_[ni].size());j++ ) {
            int n = ni;
            int m = nm_[ni][j];
						std::vector<std::vector<double>> rad = rad_poly_store_[counter];
            counter++;
//            std::cout << "n, m " << n << ", " << m << " " << counter << std::endl;
//            rad = rad_poly(rho,n,m);

            std::complex<double> Zsum(0);
            for (int xi=0; xi<dim_; xi++) {
                for (int yi = 0; yi < dim_; yi++) {
		                std::complex<double> imval(0, -theta_[xi][yi] * double(m) );
//                    Zprod[xi][yi] = grid2[xi][yi] * rad[xi][yi] * exp(imval);
                    Zprod[xi][yi] = grid[xi][yi] * rad[xi][yi] * exp(imval);
                    Zsum += Zprod[xi][yi];
//                  if (grid[xi][yi] > 0 )
//                       std::cout << xi << " " << yi << " " << grid[xi][yi] << " " << rad[xi][yi] << " " << theta_[xi][yi] << " " << Zprod[xi][yi] <<  std::endl;

                }
            }
            double norm_factor = double(n+1)/(count_nonzero(rho_) + 1 );
            Zsum=Zsum*norm_factor;
            double Znorm = std::norm(Zsum);
            invariants_.push_back(std::sqrt(Znorm));
            moments_.push_back(Zsum);
//            std::cout << "n, m " << n << ", " << m << " " << Zsum << std::endl;
//            std::cout << "ZT: " << Zsum << std::endl;
        }
    }
//    auto finish = std::chrono::high_resolution_clock::now();
 //   std::chrono::duration<double> elapsed = finish - start;
//    std::cout << "Time per transform: " << elapsed.count()  << std::endl;

//		TR << "Time for transform: " << elapsed.count() << std::endl;
  }

}
}
}
}
