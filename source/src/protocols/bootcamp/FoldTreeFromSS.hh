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
#include <protocols/loops/Loop.hh>

// Utility headers

/// Project headers
#include <core/types.hh>

// C++ headers

//Auto Headers
#include <core/pack/dunbrack/DunbrackRotamer.hh>

namespace protocols {
namespace bootcamp  {

class FoldTreeFromSS {
public:
   FoldTreeFromSS( std::string const &);

   core::kinematics::FoldTree const &
   fold_tree() const;

    core::kinematics::FoldTree
    fold_tree_from_dssp_string(std::string const & ss_string);

   protocols::loops::Loop const &
   loop( core::Size index ) const;

   core::Size
   loop_for_residue( core::Size seqpos ) const;

private:
   core::kinematics::FoldTree ft_;
   utility::vector1< protocols::loops::Loop > loop_vector_;
   utility::vector1< core::Size > loop_for_residue_;
};

}
}

