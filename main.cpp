/*
    $Id: main.cpp,v 1.4 2004/05/21 05:50:13 taku-ku Exp $;
 
   Copyright (C) 2004 Taku Kudo, All rights reserved.
     This is free software with ABSOLUTELY NO WARRANTY.
  
   This program is free software; you can redistribute it and/or modify
     it under the terms of the GNU General Public License as published by
     the Free Software Foundation; either version 2 of the License, or
     (at your option) any later version.
    
   This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU General Public License for more details.
    
   You should have received a copy of the GNU General Public License
     along with this program; if not, write to the Free Software
     Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
     02111-1307, USA
*/
#include "gspan.h"
#include "testability.hpp"

#include <chrono>
#include <functional>
#include <vector>
#include <utility>

#include <unistd.h>
#include <fstream>

#define OPT " [-m minsup] [-d] [-e] [-w] "

int main (int argc, char **argv)
{
  unsigned int maxpat = 0xffffffff;
  unsigned int minnodes = 0;
  bool where = false;
  bool enc = false;
  bool directed = false;

  std::string input_filename = "NCI1_10";

  std::ifstream input_f(input_filename);
  //  std::string input_labels("MUTAG_label");
  std::ofstream output_f("output");

  GSPAN::gSpan gspan(input_f, output_f, maxpat, minnodes, enc, where, directed);
  // auto res = gspan.run(minsup);
	 
  namespace th = thesis;
	 
  using namespace std::placeholders;
  using namespace ranges;

  auto alpha =  .05;
  
  auto to_run = [](auto f, auto gspan){return gspan.run(f);};
  auto to_run_m = [=](auto f, auto pv, auto gspan){return gspan.c_run_m(f,pv,alpha);};
  auto alg = [=](auto f){ return to_run(f,gspan);};
  auto alg_m = [=](auto f){return alg(f).size();};
  auto c_alg_m = [=](auto f, auto pv){return to_run_m(f,pv,gspan);};

  auto n1 = 15u;
  auto n2 = 44u;


  auto start = std::chrono::system_clock::now();  

  auto out = th::one_pass_(alg, n1, n2, alpha);
  //  std::cout << "alpha: " << alpha << std::endl;
  //auto out = th::bis_leap_(c_alg_m, n1, n2, alpha);

  //  auto app = th::one_pass_(alg,n1, n2, alpha);
  //  std::cout << alg_m(14) << std::endl;
  std::cout << "ROOT FREQ: " << out << std::endl;
  out = th::early_term_(c_alg_m,n1, n2, alpha);

  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> diff = end - start;
  auto test_no = alg_m(out);

  //  std::cout << "****************************************\n" << input_filename << std::endl;
  std::cout << "ROOT FREQ: " << out << std::endl;
  // std::cout << "#TESTABLES: " << test_no << std::endl;
  // std::cout << "DELTA TAR: " << alpha / test_no << std::endl;
  // std::cout << "TIME: " << diff.count() << std::endl;
  // std::cout << "DELTA BONF: " << alpha / alg_m(1) << std::endl;

  //std::cout << alg_m(14) << std::endl;

  // std::cout << "CHECKING\n";
  // std::cout << alg_m(out) << " | "  << alg_m(out+1) << " | " << alg_m(out+2) << std::endl
    ;

  //  std::cout << "alpha: " << alpha << std::endl;
  
  //  std::cout << "APP: " << app << std::endl;



}
