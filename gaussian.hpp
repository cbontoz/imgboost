/**
 * \file gaussian.hpp
 * \brief Part of imgboost library, to support gaussian/blur image filtering.
 * It offers options for (1)grayscale or RGB output and (2)'naive'
 * implementation using 2D kernel or the equivalent separated kernel in two
 * vectors. Also, the user can apply different standard deviation on X and Y
 * direction and rotate the kernel (in rad).
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

#ifndef _GAUSSIAN_HPP_
#define _GAUSSIAN_HPP_

class gaussian_abstract {
public:
  imgboost perform(imgboost const &img, size_t KernelSize, double Std);
  imgboost perform(imgboost const &img, size_t KernelSize, double StdX,
                   double StdY);

protected:
  static constexpr double _pi = M_PI;
  virtual imgboost perform(imgboost const &img, size_t KernelSize, double Rad,
                           double StdX, double StdY) = 0;
};

class gaussian : public gaussian_abstract {
public:
  imgboost perform(imgboost const &img, size_t KernelSize, double Rad,
                   double StdX, double StdY) override;

private:
  smartix<double> Create2DK(size_t KernelSize, double Rad, double StdX,
                            double StdY);
};

class gaussian_fast : public gaussian_abstract {
public:
  imgboost perform(imgboost const &img, size_t KernelSize, double Rad,
                   double StdX, double StdY) override;

private:
  smartix<double> CreateVK(size_t KernelSize, double Std);
};

imgboost gaussian_abstract::perform(imgboost const &img, size_t KernelSize,
                                    double Std) {
  return perform(img, KernelSize, 0, Std, Std);
}

imgboost gaussian_abstract::perform(imgboost const &img, size_t KernelSize,
                                    double StdX, double StdY) {
  return perform(img, KernelSize, 0, StdX, StdY);
}

imgboost gaussian::perform(imgboost const &img, size_t KernelSize, double Rad,
                           double StdX, double StdY) {

  if (KernelSize % 2 == 0)
    KernelSize += 1;

  size_t ks = KernelSize;
  size_t hks = floor(ks / 2);
  smartix<double> fw = Create2DK(KernelSize, Rad, StdX, StdY);

  imgboost ctemp = img;
  ctemp.MirrorPadding(hks);
  imgboost result = ctemp;

  size_t ctemph = ctemp.height(), ctempw = ctemp.width();

  double rsum;
  double gsum;
  double bsum;
  size_t w_inx, h_inx, f_inx;
  for (size_t h = hks; h < ctemph - hks; ++h)
    for (size_t w = hks; w < ctempw - hks; ++w) {
      rsum = 0;
      gsum = 0;
      bsum = 0;

      for (size_t h2 = 0; h2 < ks; ++h2)
        for (size_t w2 = 0; w2 < ks; ++w2) {
          w_inx = w - hks + w2;
          h_inx = h - hks + h2;
          f_inx = h2 * KernelSize + w2;
          rsum += ctemp.getRed(w_inx, h_inx) * fw.uacc(f_inx);
          gsum += ctemp.getGreen(w_inx, h_inx) * fw.uacc(f_inx);
          bsum += ctemp.getBlue(w_inx, h_inx) * fw.uacc(f_inx);
        }
      result.setPixel(w, h, imgboost::Pixel(rsum, gsum, bsum));
    }
  result.Crop(hks, hks);
  return result;
}

imgboost gaussian_fast::perform(imgboost const &img, size_t KernelSize,
                                double Rad, double StdX, double StdY) {
  if (KernelSize % 2 == 0)
    KernelSize += 1;

  size_t ks = KernelSize;
  ssize_t hks = floor(ks / 2);
  smartix<double> fwX = CreateVK(ks, StdX);
  smartix<double> fwY = CreateVK(ks, StdY);

  imgboost ctemp = img;
  ctemp.MirrorPadding(2 * hks);
  imgboost result = ctemp;
  double rsum = 0, gsum = 0, bsum = 0;
  size_t ctemph = ctemp.height(), ctempw = ctemp.width();

  std::vector<size_t> v = {ctemp.height(), ctemp.width()};
  smartix<double> cgray(v);
  smartix<double> r(v);
  smartix<double> g(v);
  smartix<double> b(v);
  for (size_t h = 0; h < ctemph; ++h)
    for (size_t w = 0; w < ctempw; ++w) {
      cgray.uacc(ctemp.findcell(h, w)) = ctemp.getGray(w, h);
      r.uacc(ctemp.findcell(h, w)) = ctemp.getRed(w, h);
      g.uacc(ctemp.findcell(h, w)) = ctemp.getGreen(w, h);
      b.uacc(ctemp.findcell(h, w)) = ctemp.getBlue(w, h);
    }

  int xhoriz[ks], xvert[ks], yhoriz[ks], yvert[ks];
  for (size_t i = 0; i < ks; ++i) {
    xhoriz[i] = (floor(((float)i - (float)hks) * cos(Rad) + 0.5));
    xvert[i] = (floor(-((float)i - (float)hks) * sin(Rad) + 0.5));
    yhoriz[i] = (floor(-((float)i - (float)hks) * sin(Rad) + 0.5));
    yvert[i] = (floor(-((float)i - (float)hks) * cos(Rad) + 0.5));
  }

  for (size_t h = hks; h < ctemph - hks; ++h)
    for (size_t w = 2 * hks; w < ctempw - 2 * hks; ++w) {
      rsum = 0;
      gsum = 0;
      bsum = 0;
      for (size_t i = 0; i < ks; ++i) {
        rsum +=
            (double)ctemp.getRed(w + xhoriz[i], h + yhoriz[i]) * fwX.uacc(i);
        gsum +=
            (double)ctemp.getGreen(w + xhoriz[i], h + yhoriz[i]) * fwX.uacc(i);
        bsum +=
            (double)ctemp.getBlue(w + xhoriz[i], h + yhoriz[i]) * fwX.uacc(i);
      }
      r.uacc(ctemp.findcell(h, w)) = rsum;
      g.uacc(ctemp.findcell(h, w)) = gsum;
      b.uacc(ctemp.findcell(h, w)) = bsum;
    }

  for (size_t h = 2 * hks; h < ctemph - 2 * hks; ++h)
    for (size_t w = hks; w < ctempw - hks; ++w) {
      rsum = 0;
      gsum = 0;
      bsum = 0;
      for (size_t i = 0; i < ks; ++i) {
        rsum +=
            r.uacc(ctemp.findcell(h + yvert[i], w + xvert[i])) * fwY.uacc(i);
        gsum +=
            g.uacc(ctemp.findcell(h + yvert[i], w + xvert[i])) * fwY.uacc(i);
        bsum +=
            b.uacc(ctemp.findcell(h + yvert[i], w + xvert[i])) * fwY.uacc(i);
      }
      result.setPixel(w, h, imgboost::Pixel(rsum, gsum, bsum));
    }
  result.Crop(2 * hks, 2 * hks);
  return result;
}

smartix<double> gaussian::Create2DK(size_t KernelSize, double Rad, double StdX,
                                    double StdY) {
  size_t ks = KernelSize;
  std::vector<size_t> v = {ks, ks};
  smartix<double> fw(v);

  double xth, yth, fsum;
  ssize_t ylm = floor(ks / 2);
  ssize_t xlm = floor(ks / 2);

  fsum = 0;
  const double _gX = 1.0 / sqrt(2.0 * _pi * pow(StdX, 2));
  const double _devX = 2.0 * pow(StdX, 2);
  const double _gY = 1.0 / sqrt(2.0 * _pi * pow(StdY, 2));
  const double _devY = 2.0 * pow(StdY, 2);
  for (ssize_t y = -ylm; y <= ylm; y++) {
    for (ssize_t x = -xlm; x <= xlm; x++) {
      xth = x * cos(Rad) - y * sin(Rad);
      yth = y * cos(Rad) + x * sin(Rad);
      fw.uacc((y + ylm) * ks + x + xlm) =
          _gY * _gX * exp(-(pow(xth, 2) / _devX + pow(yth, 2) / _devY));
      fsum += fw.uacc((y + ylm) * ks + x + xlm);
    }
  }

  // Normalise kernel
  for (ssize_t y = -ylm; y <= ylm; y++) {
    for (ssize_t x = -xlm; x <= xlm; x++) {
      fw.uacc((y + ylm) * ks + x + xlm) /= fsum;
    }
  }
  return fw;
}

smartix<double> gaussian_fast::CreateVK(size_t KernelSize, double Std) {
  size_t ks = KernelSize;
  ssize_t hks = floor(ks / 2);
  std::vector<size_t> v = {ks};
  smartix<double> fw(v);
  const double _g = 1.0 / sqrt(2.0 * _pi * pow(Std, 2));
  const double _dev = 2.0 * pow(Std, 2);
  double sum = 0;

  for (ssize_t i = -hks; i <= hks; ++i) {
    fw.uacc(i + hks) = _g * exp(-pow(i, 2) / _dev);
    sum += fw.uacc(i + hks);
  }

  for (ssize_t i = -hks; i <= hks; ++i) {
    fw.uacc(i + hks) /= sum;
  }

  return fw;
}

#endif //_GAUSSIAN_HPP_
