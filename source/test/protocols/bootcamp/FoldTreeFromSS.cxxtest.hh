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

#include <utility/vector1.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/dssp/Dssp.hh>

// Utility headers

/// Project headers
#include <core/types.hh>

// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>


using namespace protocols::match;
using namespace protocols::match::upstream;


// --------------- Test Class --------------- //

class ProteinSCSamplerTests : public CxxTest::TestSuite {

public:


	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.


	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	// --------------- Test Cases --------------- //
    void test_hello_world() {
    TS_ASSERT( true );
    }

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
    
    void test_SS_Span1() {
        utility::vector1< char > str1 = "   EEEEE   HHHHHHHH  EEEEE   IGNOR EEEEEE   HHHHHHHHHHH  EEEEE  HHHH   ";
        utility::vector1< std::pair< core::Size, core::Size > > output1 =
            identify_secondary_structure_spans(str1);
        
        utility::vector1 <core::Size> simpleArray = 0;
        utility::vector1 <core::Size> correctAnswer = {4, 8, 12, 19, 22, 26, 36, 41, 45, 55, 58, 62, 65, 68};
        
        for (core::Size ele = 1, ele < output1.size() + 1, ++ele) {
            simpleArray.push_back(outpu1[ele].first);
            simpleArray.push_back(outpu1[ele].second);
        }
        
        
        TS_ASSERT_EQUAL(simpleArray, correctAnswer);
    }

};


//Foldtree:
core::kinematics::FoldTree
fold_tree_from_ss(core::pose::PoseOP const &mypose) {
    core::scoring::dssp::Dssp dssp1(mypose);
    std::string str1 = dssp1.get_dssp_secstruct();
    return fold_tree_from_dssp_string(str1);
}


core::kinematics::FoldTree
fold_tree_from_dssp_string(std::string const & ss_string) {
    utility::vector1< std::pair< core::Size, core::Size > > output1 =
        identify_secondary_structure_spans(ss_string);
    core::kinematics::FoldTree ft;
    core::Size totE = output1.size() * 6 - 4;
    
    utility::vector1 <core::Size> segCenters = 0;
    utility::vector1 <core::Size> loopCenters = 0;
    utility::vector1 <core::Size> segbound = output1;
    
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
    
    return ft;
    
}
