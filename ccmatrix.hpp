/**
 * \file ccmatrix.hpp
 * \brief Part of imgboost library, to support co-occurrence matrix calculation.
 *
 * MIT License
 *
 * Copyright (c) 2017-2018 Christos Bontozoglou
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * \author Christos Bontozoglou
 * \email cbontoz@gmail.com
 */

#include <cassert>
#include <fstream>
#include <imgboost.hpp>
#include <iostream>
#include <math.h>
#include <memory>
#include <smartix.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <tiffio.h>
#include <vector>

#ifndef _CCMATRIX_HPP_
#define _CCMATRIX_HPP_

class ccmatrix {
private:
  size_t _lvls, _dang, _disp;
  double _amax, _amin, _astep;
  size_t _i_amax;
  double _dmin, _dmax, _dstep;
  size_t _i_dmax;
  smartix<int32_t> _ofs;
  int32_t &_offset(size_t xy01, size_t d, size_t th) {
    return _ofs.uacc(th * _i_dmax * 2 + d * 2 + xy01);
  }
  smartix<double> _ccm;
  double &_ccmatrix(size_t d, size_t a, size_t i, size_t j) {
    return _ccm.uacc(a * _i_dmax * _lvls * _lvls + d * _lvls * _lvls +
                     i * _lvls + j);
  }
  smartix<double> _fenergy;
  smartix<double> _fentropy;
  smartix<double> _fcontrast;
  smartix<double> _fhomogeneity;
  smartix<double> _fcorrelation;
  double &_feature(smartix<double> &f, size_t d, size_t a) {
    if (d >= _i_dmax)
      throw std::invalid_argument("Requested d greater than max displacement.");
    if (a >= _i_amax)
      throw std::invalid_argument("Requested a greater than max angle.");
    return f.uacc(a * _i_dmax + d);
  }

public:
  ccmatrix(size_t levels, size_t angular_resolution_degrees,
           size_t displacement);
  smartix<double> perform(imgboost &img, std::vector<double> levels);
  smartix<double> perform(imgboost &img, double low_limit, double high_limit);
  double getEnergy(size_t displacement, size_t rad);
  double getEntropy(size_t displacement, size_t rad);
  double getContrast(size_t displacement, size_t rad);
  double getHomogeneity(size_t displacement, size_t rad);
  double getCorrelation(size_t displacement, size_t rad);

  void printlast(void);
  void printfeatures(void);
  void printavgfeatures(void);
};

ccmatrix::ccmatrix(size_t levels, size_t angular_steps, size_t displacement)
    : _lvls(levels), _dang(angular_steps), _disp(displacement), _amax(M_PI),
      _amin(0), _astep((_amax - _amin) / _dang),
      _i_amax(round((_amax - _amin) / _astep)), _dmin(1), _dmax(_disp),
      _dstep(1), _i_dmax(round((_dmax - _dmin) / _dstep)),
      _ofs(2, _i_dmax, _i_amax), _ccm(_i_dmax, _i_amax, _lvls, _lvls),
      _fenergy(_i_dmax, _i_amax), _fentropy(_i_dmax, _i_amax),
      _fcontrast(_i_dmax, _i_amax), _fhomogeneity(_i_dmax, _i_amax),
      _fcorrelation(_i_dmax, _i_amax) {
  (_dang < 1) ? throw std::invalid_argument("Angular Resolution must be > 0")
              : 0;
  (_disp < 1) ? throw std::invalid_argument("Displacement should be > 0") : 0;

  // Calculate displacements' tables
  for (size_t thi = 0; thi < _i_amax; ++thi) {
    for (size_t di = 0; di < _i_dmax; ++di) {
      _offset(0, di, thi) = static_cast<double>(di + _dmin) * _dstep *
                            sin(static_cast<double>(thi) * _astep);
      _offset(1, di, thi) = static_cast<double>(di + _dmin) * _dstep *
                            cos(static_cast<double>(thi) * _astep);
    }
  }
}

smartix<double> ccmatrix::perform(imgboost &img, std::vector<double> limits) {
  smartix<uint8_t> limg(img.scalar());
  if ((_lvls + 1) != limits.capacity())
    throw std::invalid_argument("CMatrix levels do not match levels' vector");
  // Replace array values with levels
  for (size_t h = 0; h < img.height(); ++h)
    for (size_t w = 0; w < img.width(); ++w)
      for (size_t l = 0; l < limits.size() - 1; ++l) {
        if (img.scalar(w, h) >= limits[l] && img.scalar(w, h) < limits[l + 1])
          limg.uacc(h * img.width() + w) = l;
        else if (img.scalar(w, h) >= limits[limits.size() - 1])
          limg.uacc(h * img.width() + w) = l - 1;
      }

  // calc matrix
  _ccm = 0;
  _fenergy = 0;
  _fentropy = 0;
  _fcontrast = 0;
  _fhomogeneity = 0;
  _fcorrelation = 0;
  for (size_t radi = 0; radi < _i_amax; ++radi) {
    for (size_t di = 0; di < _i_dmax; ++di) {
      int32_t offx, offy;
      offx = _offset(1, di, radi);
      offy = _offset(0, di, radi);

      int32_t ystr, ystp;
      if (offy < 0) {
        ystr = -offy;
        ystp = img.height();
      } else {
        ystr = 0;
        ystp = img.height() - offy;
      }
      int32_t xstr, xstp;
      if (offx < 0) {
        xstr = -offx;
        xstp = img.width();
      } else {
        xstr = 0;
        xstp = img.width() - offx;
      }

      for (int y = ystr; y < ystp; y++) {
        for (int x = xstr; x < xstp; x++) {
          ++_ccmatrix(di, radi, limg.uacc(y * _lvls + x),
                      limg.uacc((y + offy) * _lvls + (x + offx)));
          ++_ccmatrix(di, radi, limg.uacc((y + offy) * _lvls + (x + offx)),
                      limg.uacc(y * _lvls + x));
        }
      }

      double sum = 0;
      for (size_t i = 0; i < _lvls; ++i)
        for (size_t j = 0; j < _lvls; ++j) {
          sum += _ccmatrix(di, radi, i, j);
        }

      double mean = 0, imj2;
      for (size_t i = 0; i < _lvls; ++i)
        for (size_t j = 0; j < _lvls; ++j) {
          _ccmatrix(di, radi, i, j) /= sum;
          _feature(_fenergy, di, radi) +=
              _ccmatrix(di, radi, i, j) * _ccmatrix(di, radi, i, j);
          if (_ccmatrix(di, radi, i, j) != 0.0)
            _feature(_fentropy, di, radi) +=
                (-log(_ccmatrix(di, radi, i, j))) * _ccmatrix(di, radi, i, j);
          imj2 = pow(static_cast<double>(i) - static_cast<double>(j), 2);
          _feature(_fcontrast, di, radi) += _ccmatrix(di, radi, i, j) * imj2;
          _feature(_fhomogeneity, di, radi) +=
              _ccmatrix(di, radi, i, j) / (1.0 + imj2);
          mean += static_cast<double>(i) * _ccmatrix(di, radi, i, j);
        }
      double std = 0;
      for (size_t i = 0; i < _lvls; ++i)
        for (size_t j = 0; j < _lvls; ++j) {
          std +=
              _ccmatrix(di, radi, i, j) * pow(static_cast<double>(i) - mean, 2);
        }
      for (size_t i = 0; i < _lvls; ++i)
        for (size_t j = 0; j < _lvls; ++j) {
          _feature(_fcorrelation, di, radi) +=
              _ccmatrix(di, radi, i, j) * (static_cast<double>(i) - mean) *
              (static_cast<double>(j) - mean) / pow(std, 2);
        }
    }
  }
  return _ccm;
}

smartix<double> ccmatrix::perform(imgboost &img, double low_limit,
                                  double high_limit) {
  std::vector<double> limits;
  double lstep = (high_limit - low_limit) / _lvls;
  for (size_t l = 0; l < _lvls + 1; ++l) {
    limits.push_back(static_cast<double>(l) * lstep + low_limit);
  }
  return perform(img, limits);
}

void ccmatrix::printlast(void) {
  std::ofstream myfile;
  myfile.open("./example.csv");
  for (size_t radi = 0; radi < _i_amax; ++radi)
    for (size_t di = 0; di < _i_dmax; ++di) {
      myfile << "Angle:," << static_cast<double>(radi) * _astep * 180.0 / M_PI
             << ", Displacement:," << di + _dmin << "\n";
      for (size_t i = 0; i < _lvls; ++i) {
        for (size_t j = 0; j < _lvls; ++j) {
          myfile << _ccmatrix(di, radi, i, j) << ",";
        }
        myfile << "\n";
      }
      myfile << "\n";
    }
  myfile.close();
}

void ccmatrix::printfeatures(void) {
  std::ofstream myfile;
  myfile.open("./features.csv");
  for (size_t radi = 0; radi < _i_amax; ++radi)
    for (size_t di = 0; di < _i_dmax; ++di) {
      myfile << "Angle:," << static_cast<double>(radi) * _astep * 180.0 / M_PI
             << ",Displacement:," << di + _dmin << "\n";
      myfile << "Energy,Entropy,Contrast,Homogeneity,Correlation\n";
      myfile << _feature(_fenergy, di, radi) << ","
             << _feature(_fentropy, di, radi) << ","
             << _feature(_fcontrast, di, radi) << ","
             << _feature(_fhomogeneity, di, radi) << ","
             << _feature(_fcorrelation, di, radi) << "\n";
      myfile << "\n";
    }
  myfile.close();
}

void ccmatrix::printavgfeatures(void) {
  std::ofstream myfile;
  myfile.open("./avgfeatures.csv");
  myfile << "Rad,Degrees,Energy,Entropy,Contrast,Homogeneity,Correlation\n";
  for (size_t radi = 0; radi < _i_amax; ++radi) {
    double seng = 0, sent = 0, scon = 0, shomo = 0, scor = 0;
    for (size_t di = 0; di < _i_dmax; ++di) {
      seng += _feature(_fenergy, di, radi);
      sent += _feature(_fentropy, di, radi);
      scon += _feature(_fcontrast, di, radi);
      shomo += _feature(_fhomogeneity, di, radi);
      scor += _feature(_fcorrelation, di, radi);
    }
    myfile << static_cast<double>(radi) * _astep << ","
           << static_cast<double>(radi) * _astep * 180.0 / M_PI << ","
           << seng / static_cast<double>(_i_dmax) << ","
           << sent / static_cast<double>(_i_dmax) << ","
           << scon / static_cast<double>(_i_dmax) << ","
           << shomo / static_cast<double>(_i_dmax) << ","
           << scor / static_cast<double>(_i_dmax) << "\n";
  }
  myfile.close();
}

double ccmatrix::getEnergy(size_t displacement, size_t rad) {
  size_t a_inx = round((rad - _amin) / _astep);
  return _feature(_fenergy, displacement, a_inx);
}
double ccmatrix::getEntropy(size_t displacement, size_t rad) {
  size_t a_inx = round((rad - _amin) / _astep);
  return _feature(_fentropy, displacement, a_inx);
}
double ccmatrix::getContrast(size_t displacement, size_t rad) {
  size_t a_inx = round((rad - _amin) / _astep);
  return _feature(_fcontrast, displacement, a_inx);
}

double ccmatrix::getHomogeneity(size_t displacement, size_t rad) {
  size_t a_inx = round((rad - _amin) / _astep);
  return _feature(_fhomogeneity, displacement, a_inx);
}

double ccmatrix::getCorrelation(size_t displacement, size_t rad) {
  size_t a_inx = round((rad - _amin) / _astep);
  return _feature(_fcorrelation, displacement, a_inx);
}

#endif //_CCMATRIX_HPP_
