/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2012-2014  Universite catholique de Louvain (UCL), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <utility>
#include <vector>

#include "TApplication.h"
#include "TROOT.h"
#include "TSystem.h"

#include "TString.h"

#include "TClonesArray.h"
#include "TH2.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TPaveText.h"

#include "classes/DelphesClasses.h"

#include "ExRootAnalysis/ExRootResult.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeReader.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#include "ExRootAnalysis/ExRootUtilities.h"

using namespace std;

//------------------------------------------------------------------------------

// Here you can put your analysis macro

#include "Example1.C"

//------------------------------------------------------------------------------

void usage(){
            std::cerr << " Usage: Example1 <option(s)> SOURCES"
              << "Options:\n"
              << "\t-h,--help\t\tShow this help message\n"
              << "\t-i,--input\tSpecify the input filename\n"
              << "\t-o,--output\tSpecify the output filename\n"
              << "\t-d,--debug\tEnable debug\n"
              << "\t-g,--genJet\tEnable GenJet\n"
              << std::endl;
}
int main(int argc, char *argv[])
{
  char *appName = "Example1";


  TString inputFile(""), outputFile("");
  bool isDebug=false;
  bool isGenJet=false;

  for (int i = 1; i < argc; i++) {
        
        std::string arg = argv[i];

        if ((arg == "-h") || (arg == "--help")) {
            usage();
            return 0;
        } else if ((arg == "-i") || (arg == "--input")) {
            i++;
            if (i < argc) inputFile = TString(argv[i]);
            else {
             std::cerr << "--input option requires one argument." << std::endl;
                return 1;

            }
        }else if ((arg == "-o") || (arg == "--output")) {
            i++;
            if (i < argc) outputFile = TString(argv[i]);
            else {
             std::cerr << "--output option requires one argument." << std::endl;
                return 1;

            }
        } else if ((arg == "-d") || (arg == "--debug")) {
            isDebug = true;
        } else if ((arg == "-g") || (arg == "--genJet")) {
            isGenJet = true;
        }
  }

  if(inputFile==TString("") || outputFile==TString("") ){
    usage();
    return 1;
  } 

  gROOT->SetBatch();

  int appargc = 1;
  char *appargv[] = {appName};
  TApplication app(appName, &appargc, appargv);


  //------------------------------------------------------------------------------

  // Here you call your macro's main function

  Example1(inputFile.Data(), outputFile.Data(), isGenJet, isDebug);

  //------------------------------------------------------------------------------
}
