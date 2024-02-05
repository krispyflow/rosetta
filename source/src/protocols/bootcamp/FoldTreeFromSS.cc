// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/protocols/match/ProteinSCSampler.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/match/upstream/ProteinSCSampler.hh>
#include <protocols/match/upstream/OriginalScaffoldBuildPoint.hh>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>


#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <protocols/loops/Loop.hh>

// Utility headers
#include <utility/vector1.hh>
/// Project headers
#include <core/types.hh>

// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <protocols/bootcamp/FoldTreeFromSS.hh>

namespace protocols {
namespace bootcamp  {

FoldTreeFromSS::
FoldTreeFromSS( core::pose::PoseOP const &mypose ) {
  core::scoring::dssp::Dssp dssp1(mypose);
  std::string str1 = dssp1.get_dssp_secstruct();
  FoldTreeFrom(str1);
}

FoldTreeFromSS::
FoldTreeFromSS( std::string const & ss_string ) {
   utility::vector1< std::pair< core::Size, core::Size > > output1 =
        identify_secondary_structure_spans(ss_string);
    core::kinematics::FoldTree ft;
    core::Size totE = output1.size() * 6 - 4;
    
    utility::vector1 <core::Size> segCenters;
    utility::vector1 <core::Size> loopCenters;
    utility::vector1< std::pair< core::Size, core::Size > >  segbound = output1;
    
    segbound[1].first = 1;
    segbound[-1].second = ss_string.size();
    
    for (core::Size E = 1; E <= output1.size(); ++E) {
        //find 2D segment nodes.
        segCenters.push_back((output1[E].first + output1[E].second) / 2);
        if (E != totE) {
            loopCenters.push_back((output1[E].second + output1[E + 1].first)/2);
        }
    }
    
    segbound[1] = 1;
    segbound[-1] = ss_string.size();
    
    //build jump edgies:
    core::Size totJump = output1.size() * 2 - 2;
    for (core::Size jE = 1; jE <= totJump; ++jE) {
        ft.add_edge(segCenters[1], loopCenters[jE], jE);
    }
    
    //build 2D structure peptide edgies:
    for (core::Size pE = 1; pE <= output1.size(); ++pE) {
        ft.add_edge(segCenter[pE], segbound[pE].first, core::kinematics::Edge::PEPTIDE);
        ft.add_edge(segCenter[pE], segbound[pE].second, core::kinematics::Edge::PEPTIDE);
    }
    
    //build loop peptide edgies:
    for (core::Size pE = 1; pE <= output1.size() - 1; ++pE) {
        ft.add_edge(loopCenter[pE], segbound[pE].second + 1, core::kinematics::Edge::PEPTIDE);
        ft.add_edge(loopCenter[pE], segbound[pE + 1].first - 1, core::kinematics::Edge::PEPTIDE);
    }
    ft_ = ft;
    loop_for_residue_ = segbound;

    //make a vector of nodes posi only, should not need sorting
    utility::vector1 <core::Size> allCenters;
    for (core::Size i = 1; i <= segCenters.size(); ++i) {
      allCenters.push_back(segCenters[i]);
      if (i != segCenters.size()) {
        allCenters.push_back(loopCenters[i]);
      }
    }  

    //make a vector for loop_for_residue_
    utility::vector1 <core::Size> loop_for_residue;
    for (core::Size resiItr = 1; resiItr <= ss_string.size(); ++resiItr) {
      core::Size index = 0;
      core::Size nodeIndex = 1;
      if (resiItr < allCenters[nodeIndex]) {
        loop_for_residue.push_back(index);
      } else if (resiItr < allCenters[-1]) {
        ++index;
        ++nodeIndex;
        loop_for_residue.push_back(index);
      } else {
        loop_for_residue.push_back(0);
      }
    }
    loop_for_residue_ = loop_for_residue;

    //how many loops we have (== num of jump edgies)??
    //first make a vector of boundary residue:
    utility::vector1 <core::Size> bound;
    for (core::Size i = 1; i <= segbound.size(); ++i) {
      bound.push_back(segbound[i].first);
      bound.push_back(segbound[i].second);
    }
    bound.erase(bound.begin());
    bound.pop_back();

    utility::vector1< protocols::loops::Loop > loop_vector;
    for (core::Size i = 1; i <= totJump; ++i) {
      loop_vector.push_back(protocols::loop::Loop(allCenters[i], allCenters[i+1], bound[i]));
    } 
    loop_vector_ = loop_vector;
}

core::kinematics::FoldTree const &
fold_tree() const{
  return ft_;
}

protocols::loops::Loop const &
loop( core::Size index ) const{
  return loop_vector_.[index];
}

//returns cut-point
core::Size
loop_for_residue( core::Size seqpos ) const{
  return loop_for_residue_.[seqos];
}

// private:
//    core::kinematics::FoldTree ft_;
//    utility::vector1< protocols::loops::Loop > loop_vector_;
//    utility::vector1< core::Size > loop_for_residue_;
// };

utility::vector1< std::pair< core::Size, core::Size > >
identify_secondary_structure_spans( std::string const & ss_string )
{
  utility::vector1< std::pair< core::Size, core::Size > > ss_boundaries;
  core::Size strand_start = -1;
  for ( core::Size ii = 0; ii < ss_string.size(); ++ii ) {
    if ( ss_string[ ii ] == 'E' || ss_string[ ii ] == 'H'  ) {
      if ( int( strand_start ) == -1 ) {
        strand_start = ii;
      } else if ( ss_string[ii] != ss_string[strand_start] ) {
        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
        strand_start = ii;
      }
    } else {
      if ( int( strand_start ) != -1 ) {
        ss_boundaries.push_back( std::make_pair( strand_start+1, ii ) );
        strand_start = -1;
      }
    }
  }
  if ( int( strand_start ) != -1 ) {
    // last residue was part of a ss-eleemnt
    ss_boundaries.push_back( std::make_pair( strand_start+1, ss_string.size() ));
  }
  for ( core::Size ii = 1; ii <= ss_boundaries.size(); ++ii ) {
    std::cout << "SS Element " << ii << " from residue "
      << ss_boundaries[ ii ].first << " to "
      << ss_boundaries[ ii ].second << std::endl;
  }
  return ss_boundaries;
}


//Foldtree:
// core::kinematics::FoldTree
// fold_tree_from_ss(core::pose::PoseOP const &mypose) {
//     core::scoring::dssp::Dssp dssp1(mypose);
//     std::string str1 = dssp1.get_dssp_secstruct();
//     return fold_tree_from_dssp_string(str1);
// }


// core::kinematics::FoldTree
// fold_tree_from_dssp_string(std::string const & ss_string) {
//     utility::vector1< std::pair< core::Size, core::Size > > output1 =
//         identify_secondary_structure_spans(ss_string);
//     core::kinematics::FoldTree ft;
//     core::Size totE = output1.size() * 6 - 4;
    
//     utility::vector1 <core::Size> segCenters = 0;
//     utility::vector1 <core::Size> loopCenters = 0;
//     utility::vector1 <core::Size> segbound = output1;
    
//     segbound[1].first = 1;
//     segbound[-1].second = ss_string.size();
    
//     for (core::Size E = 1; E <= output1.size(); ++E) {
//         //find 2D segment nodes.
//         segCenters.push_back((output1[E].first + output1[E].second) / 2);
//         if (E != totE) {
//             loopCenters.push_back((output1[E].second + output1[E + 1].first)/2);
//         }
//     }
    
//     segbound[1] = 1;
//     segbound[-1] = ss_string.size();
    
//     //build jump edgies:
//     core::Size totJump = output1.size() * 2 - 2;
//     for (core::Size jE = 1; jE <= totJump; ++jE) {
//         ft.add_edge(segCenters[1], loopCenters[jE], jE);
//     }
    
//     //build 2D structure peptide edgies:
//     for (core::Size pE = 1; pE <= output1.size(); ++pE) {
//         ft.add_edge(segCenter[pE], segbound[pE].first, core::kinematics::Edge::PEPTIDE);
//         ft.add_edge(segCenter[pE], segbound[pE].second, core::kinematics::Edge::PEPTIDE);
//     }
    
//     //build loop peptide edgies:
//     for (core::Size pE = 1; pE <= output1.size() - 1; ++pE) {
//         ft.add_edge(loopCenter[pE], segbound[pE].second + 1, core::kinematics::Edge::PEPTIDE);
//         ft.add_edge(loopCenter[pE], segbound[pE + 1].first - 1, core::kinematics::Edge::PEPTIDE);
//     }
//     return ft;
// }



}
}
